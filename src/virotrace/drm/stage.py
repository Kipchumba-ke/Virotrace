"""The drm-call stage: projected pol alignment -> mutation table and gene coverage.

Outputs (per scope):
    results/<scope>/drm/gene_coverage.tsv   one row per (strain, gene): how much of the
                                            gene the sequence actually covers, and why a
                                            gene was disqualified (frameshift, stop)
    results/<scope>/drm/drm_table.tsv       mutations at positions on the SDRM list
    results/<scope>/drm/drm_summary.json    counts, coverage deciles, the SDRM version

This stage reads every aligned sequence. It does not apply QC, subsampling or the clone
policy: those are decisions about which rows to *count*, made later against the
coverage table, so changing one never requires re-calling mutations.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from ..align.project import Insertion
from ..errors import ContractError
from ..io.fasta import read_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .caller import SequenceCalls, call_projected, insertions_by_gene
from .mutations import GENE_PROTEIN_LENGTH
from .reference import gene_proteins, load_reference_sequence
from .sdrm import load_sdrm_list

COVERAGE_COLUMNS = [
    "strain",
    "gene",
    "coverage",
    "n_called",
    "n_deleted",
    "n_too_ambiguous",
    "frameshift",
    "n_internal_stops",
    "usable",
    "identity",
    "sdrm_position_coverage",
    "sdrm_position_coverage_sensitivity",
]
TABLE_COLUMNS = [
    "strain",
    "gene",
    "position",
    "ref_aa",
    "residues",
    "is_mixture",
    "is_sdrm",
    "sdrm_residues",
    "sdrm_list_version",
    "is_sdrm_sensitivity",
    "sensitivity_residues",
]


def _decile(fraction: float) -> int:
    return min(int(fraction * 10), 9)


def nucleotide_identity(sequence: str, reference: str) -> float | None:
    """Fraction of unambiguous aligned bases that match the reference.

    The cheapest way to tell a real pol sequence (typically 80-95% identical to HXB2)
    from a sequence an aligner forced against the wrong gene (near 25-40%). ``None``
    when nothing was sequenced, so an uncovered gene is never reported as identity 0.
    """
    informative = matches = 0
    for base, ref_base in zip(sequence, reference, strict=True):
        if base in "ACGT":
            informative += 1
            matches += base == ref_base
    return matches / informative if informative else None


def run_drm_call(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    genes_cfg = settings.reference.genes
    pol_start = genes_cfg.pol[0]
    spans: dict[str, tuple[int, int]] = {
        "PR": genes_cfg.pr,
        "RT": genes_cfg.rt,
        "IN": genes_cfg.in_,
    }

    interim = root / settings.paths["interim"] / scope
    aln_path, ins_path = interim / "aln.pol.unmasked.fasta", interim / "insertions.tsv"
    for path in (aln_path, ins_path):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace align` first")

    reference_dir = root / settings.paths["reference"]
    reference_protein = gene_proteins(
        load_reference_sequence(reference_dir / settings.reference.record), spans
    )
    sdrm = load_sdrm_list(
        reference_dir / settings.drm.sdrm.list, version=settings.drm.sdrm.list_version
    )
    reference_pol = load_reference_sequence(reference_dir / settings.reference.record)[
        pol_start - 1 : genes_cfg.pol[1]
    ]
    sensitivity = (
        load_sdrm_list(
            reference_dir / settings.drm.sdrm.sensitivity_list,
            version=settings.drm.sdrm.sensitivity_list_version or "sensitivity",
        )
        if settings.drm.sdrm.sensitivity_list
        else None
    )
    sens_positions = {
        gene: {pos for (g, pos) in (sensitivity.entries if sensitivity else {}) if g == gene}
        for gene in spans
    }
    sdrm_positions = {gene: {pos for (g, pos) in sdrm.entries if g == gene} for gene in spans}
    policy = settings.drm.mixtures.policy
    max_residues = settings.drm.mixtures.max_ambiguous_residues
    min_coverage = settings.qc.min_gene_coverage

    insertions: dict[str, list[Insertion]] = defaultdict(list)
    for row in read_tsv(ins_path):
        insertions[row["strain"]].append(Insertion(int(row["after_ref_pos"]), int(row["length"])))

    coverage_rows: list[dict[str, Any]] = []
    table_rows: list[dict[str, Any]] = []
    deciles: dict[str, Counter[int]] = {gene: Counter() for gene in spans}
    meeting: Counter[str] = Counter()
    n_frameshift: Counter[str] = Counter()
    n_stop: Counter[str] = Counter()
    identity_deciles: dict[str, Counter[int]] = {gene: Counter() for gene in spans}
    covering_all: Counter[str] = Counter()
    n_sequences = n_no_pol = 0

    for strain, projected in read_fasta(interim / "aln.pol.unmasked.fasta"):
        n_sequences += 1
        result: SequenceCalls = call_projected(
            strain,
            projected,
            pol_start=pol_start,
            gene_spans=spans,
            reference_protein=reference_protein,
            insertion_lengths=insertions_by_gene(insertions.get(strain, []), spans, pol_start),
            mixture_policy=policy,
            max_ambiguous_residues=max_residues,
        )
        raw_best = 0.0
        for gene, gene_result in result.genes.items():
            length = GENE_PROTEIN_LENGTH[gene]
            raw = (len(gene_result.called) + len(gene_result.deleted)) / length
            raw_best = max(raw_best, raw)
            deciles[gene][_decile(raw)] += 1
            usable_cov = result.calls.coverage_fraction(gene)
            if usable_cov >= min_coverage:
                meeting[gene] += 1
            n_frameshift[gene] += gene_result.frameshift
            n_stop[gene] += bool(gene_result.internal_stops)
            offset = spans[gene][0] - pol_start
            identity = nucleotide_identity(
                projected[offset : offset + 3 * length], reference_pol[offset : offset + 3 * length]
            )
            if identity is not None:
                identity_deciles[gene][_decile(identity)] += 1
            wanted = sdrm_positions[gene]
            position_coverage = (
                len(wanted & result.calls.covered[gene]) / len(wanted) if wanted else None
            )
            if position_coverage == 1.0:
                covering_all[gene] += 1
            sens_wanted = sens_positions[gene]
            sens_coverage = (
                len(sens_wanted & result.calls.covered[gene]) / len(sens_wanted)
                if sens_wanted
                else None
            )
            coverage_rows.append(
                {
                    "strain": strain,
                    "gene": gene,
                    "coverage": round(usable_cov, 4),
                    "n_called": len(gene_result.called),
                    "n_deleted": len(gene_result.deleted),
                    "n_too_ambiguous": gene_result.n_too_ambiguous,
                    "frameshift": gene_result.frameshift,
                    "n_internal_stops": len(gene_result.internal_stops),
                    "usable": gene_result.usable,
                    "identity": None if identity is None else round(identity, 4),
                    "sdrm_position_coverage": (
                        None if position_coverage is None else round(position_coverage, 4)
                    ),
                    "sdrm_position_coverage_sensitivity": (
                        None if sens_coverage is None else round(sens_coverage, 4)
                    ),
                }
            )
        if raw_best < 0.10:
            n_no_pol += 1

        for mutation in result.calls.mutations:
            key = (mutation.gene, mutation.position)
            if key not in sdrm.entries and not (sensitivity and key in sensitivity.entries):
                continue
            matched = sdrm.matched_residues(mutation)
            sens_matched = sensitivity.matched_residues(mutation) if sensitivity else frozenset()
            table_rows.append(
                {
                    "strain": strain,
                    "gene": mutation.gene,
                    "position": mutation.position,
                    "ref_aa": mutation.ref_aa,
                    "residues": "".join(sorted(mutation.alt_aa)),
                    "is_mixture": mutation.is_mixture,
                    "is_sdrm": bool(matched),
                    "sdrm_residues": "".join(sorted(matched)),
                    "sdrm_list_version": sdrm.version,
                    "is_sdrm_sensitivity": bool(sens_matched),
                    "sensitivity_residues": "".join(sorted(sens_matched)),
                }
            )

    if n_sequences == 0:
        raise ContractError(f"{aln_path} contains no sequences")

    out = root / settings.paths["results"] / scope / "drm"
    write_tsv(out / "gene_coverage.tsv", coverage_rows, COVERAGE_COLUMNS)
    write_tsv(out / "drm_table.tsv", table_rows, TABLE_COLUMNS)

    summary = {
        "scope": scope,
        "n_sequences": n_sequences,
        "n_with_almost_no_pol": n_no_pol,
        "min_gene_coverage": min_coverage,
        "mixture_policy": policy,
        "sdrm_list": settings.drm.sdrm.list,
        "sdrm_list_version": sdrm.version,
        "sensitivity_list": settings.drm.sdrm.sensitivity_list,
        "n_sdrm_rows": sum(1 for r in table_rows if r["is_sdrm"]),
        "n_sequences_with_sdrm": len({r["strain"] for r in table_rows if r["is_sdrm"]}),
        "per_gene": {
            gene: {
                "n_meeting_min_coverage": meeting[gene],
                "n_frameshift": n_frameshift[gene],
                "n_internal_stop": n_stop[gene],
                "n_covering_all_sdrm_positions": covering_all[gene],
                "n_sdrm_positions_on_list": len(sdrm_positions[gene]),
                "raw_coverage_deciles": [deciles[gene][d] for d in range(10)],
                "identity_deciles": [identity_deciles[gene][d] for d in range(10)],
            }
            for gene in spans
        },
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "drm_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
