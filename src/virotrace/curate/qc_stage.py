"""The QC stage: which aligned sequences count, and why the others do not.

Runs after alignment and the per-gene coverage table, because coverage, identity,
frame and stops are all defined on HXB2 coordinates. Rules apply in order and the first
failure is the sequence's reason, so the waterfall cannot double-count.

Sequence-level rules (a failure removes the sequence):
    QC-7  laboratory strain list
    QC-0  not pol: best-gene identity to HXB2 below ``qc.min_identity_to_reference``
    QC-1  no gene is usable AND covers enough of its class's SDRM positions
    QC-2  ambiguity fraction above ``qc.max_ambiguity_frac``
    QC-5  APOBEC hypermutation (excluded or only flagged, per ``qc.hypermutation.action``)

Gene-level findings (QC-3 internal stop, QC-4 frameshift) disqualify one gene, not the
sequence: a PR+RT genotype with a bad RT is still good evidence for protease. They are
written to the ledger with the gene named.

Outputs (per scope):
    interim/<scope>/qc.passed.tsv        strains that count, with the genes that qualify
    results/<scope>/qc/excluded.tsv      the ledger: strain, gene, rule, observed, threshold
    results/<scope>/qc/qc_summary.json   the waterfall and the corpus recount
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from ..drm.reference import load_reference_sequence
from ..errors import ContractError
from ..io.fasta import read_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .clones import apply_population_policies
from .corpus import count_corpus
from .hypermut import screen_hypermutation

LEDGER_COLUMNS = ["strain", "gene", "rule", "observed", "threshold"]


def _truthy(value: str) -> bool:
    return value.strip().lower() == "true"


def _number(value: str) -> float | None:
    return float(value) if value.strip() else None


def load_lab_strains(path: Path) -> frozenset[str]:
    if not path.exists():
        return frozenset()
    lines = (line.strip() for line in path.read_text().splitlines())
    return frozenset(line for line in lines if line and not line.startswith("#"))


def qualifying_genes(
    rows: list[dict[str, str]], *, min_position_coverage: float, min_gene_coverage: float
) -> list[str]:
    """Genes that are usable and cover enough to count toward a denominator.

    Where the SDRM list has positions for the gene, the test is the share of those
    positions covered; where it has none (integrase on the 2009 list) it falls back to
    plain gene coverage.
    """
    genes = []
    for row in rows:
        if not _truthy(row["usable"]):
            continue
        position_coverage = _number(row["sdrm_position_coverage"])
        if position_coverage is not None:
            ok = position_coverage >= min_position_coverage
        else:
            ok = float(row["coverage"]) >= min_gene_coverage
        if ok:
            genes.append(row["gene"])
    return genes


def run_qc(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    qc = settings.qc
    interim = root / settings.paths["interim"] / scope
    results = root / settings.paths["results"] / scope

    coverage_path = results / "drm" / "gene_coverage.tsv"
    meta_path = interim / "metadata.valid.tsv"
    aln_path = interim / "aln.pol.unmasked.fasta"
    for path, stage in ((coverage_path, "drm call"), (meta_path, "normalise"), (aln_path, "align")):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    by_strain: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(coverage_path):
        by_strain[row["strain"]].append(row)
    metadata = {row["strain"]: row for row in read_tsv(meta_path)}
    lab = (
        load_lab_strains(root / settings.paths["reference"] / qc.lab_strain_list)
        if qc.exclude_lab_strains
        else frozenset()
    )

    ledger: list[dict[str, Any]] = []
    waterfall: Counter[str] = Counter()
    survivors: dict[str, list[str]] = {}

    def fail(strain: str, rule: str, observed: Any, threshold: Any, gene: str = "") -> None:
        ledger.append(
            {
                "strain": strain,
                "gene": gene,
                "rule": rule,
                "observed": observed,
                "threshold": threshold,
            }
        )

    for strain in by_strain:
        gene_rows = by_strain[strain]

        if strain in lab:
            fail(strain, "QC-7", "listed", "not listed")
            waterfall["QC-7"] += 1
            continue

        identities = [v for r in gene_rows if (v := _number(r["identity"])) is not None]
        best_identity = max(identities, default=None)
        if best_identity is not None and best_identity < qc.min_identity_to_reference:
            fail(strain, "QC-0", round(best_identity, 4), qc.min_identity_to_reference)
            waterfall["QC-0"] += 1
            continue

        # Gene-level findings are only meaningful for real pol: a forced mis-alignment
        # of another gene has stops and frameshifts everywhere and says nothing.
        for row in gene_rows:
            if _truthy(row["frameshift"]):
                fail(strain, "QC-4", "frameshift", "none", row["gene"])
            if int(row["n_internal_stops"]) > 0:
                fail(strain, "QC-3", row["n_internal_stops"], 0, row["gene"])

        genes = qualifying_genes(
            gene_rows,
            min_position_coverage=qc.min_class_position_coverage,
            min_gene_coverage=qc.min_gene_coverage,
        )
        if not genes:
            best = max((float(r["coverage"]) for r in gene_rows), default=0.0)
            fail(strain, "QC-1", round(best, 4), qc.min_class_position_coverage)
            waterfall["QC-1"] += 1
            continue

        ambiguity = float(metadata[strain]["ambiguity_frac"]) if strain in metadata else 0.0
        if ambiguity > qc.max_ambiguity_frac:
            fail(strain, "QC-2", round(ambiguity, 4), qc.max_ambiguity_frac)
            waterfall["QC-2"] += 1
            continue
        survivors[strain] = genes

    # QC-5 needs the sequences, so it runs only on what is left (and never on a sequence
    # that is not pol, where a context test would be meaningless).
    hypermutated: set[str] = set()
    if qc.hypermutation.enabled and survivors:
        pol_start, pol_end = settings.reference.genes.pol
        record = root / settings.paths["reference"] / settings.reference.record
        reference_pol = load_reference_sequence(record)[pol_start - 1 : pol_end]
        for strain, projected in read_fasta(aln_path):
            if strain not in survivors:
                continue
            result = screen_hypermutation(
                strain, projected, reference_pol, p_threshold=qc.hypermutation.p_threshold
            )
            if result.flagged:
                hypermutated.add(strain)
                fail(strain, "QC-5", f"p={result.p_value:.2e}", qc.hypermutation.p_threshold)
        if qc.hypermutation.action == "exclude":
            for strain in hypermutated:
                survivors.pop(strain)
            waterfall["QC-5"] += len(hypermutated)

    n_in = len(by_strain)
    passed_rows = [
        {"strain": s, "genes": ",".join(genes), "hypermutated": s in hypermutated}
        for s, genes in survivors.items()
    ]
    excluded_sequences = sum(waterfall.values())
    if len(passed_rows) + excluded_sequences != n_in:
        raise ContractError(
            f"record conservation violated in QC: {n_in} in, {len(passed_rows)} passed, "
            f"{excluded_sequences} excluded"
        )

    write_tsv(interim / "qc.passed.tsv", passed_rows, ["strain", "genes", "hypermutated"])
    write_tsv(results / "qc" / "excluded.tsv", ledger, LEDGER_COLUMNS)

    def typed(strain: str) -> dict[str, Any]:
        row: dict[str, Any] = dict(metadata[strain])
        row["year"] = int(row["year"]) if row["year"].strip() else None
        row["clone_derived"] = _truthy(row.get("clone_derived", "false"))
        return row

    kept = [typed(s) for s in survivors if s in metadata]
    target = (settings.scopes[scope].get("countries") or ["KEN"])[0]
    after_policy, removed = apply_population_policies(
        kept,
        clone_policy=qc.clone_policy,
        individual_policy=qc.individual_policy,
        patterns=qc.individual_id_patterns,
    )
    summary = {
        "scope": scope,
        "n_aligned": n_in,
        "n_passed": len(passed_rows),
        "waterfall": {
            rule: waterfall.get(rule, 0) for rule in ("QC-7", "QC-0", "QC-1", "QC-2", "QC-5")
        },
        "n_hypermutated_flagged": len(hypermutated),
        "hypermutation_action": qc.hypermutation.action,
        "thresholds": {
            "min_identity_to_reference": qc.min_identity_to_reference,
            "min_class_position_coverage": qc.min_class_position_coverage,
            "max_ambiguity_frac": qc.max_ambiguity_frac,
        },
        "lab_strain_list_entries": len(lab),
        "n_genes_qualifying": dict(Counter(g for genes in survivors.values() for g in genes)),
        "clone_policy": qc.clone_policy,
        "individual_policy": qc.individual_policy,
        "n_removed_by_population_policy": removed,
        "corpus_after_qc": count_corpus(kept, country_code=target),
        "corpus_after_qc_and_clone_policy": count_corpus(after_policy, country_code=target),
    }
    (results / "qc").mkdir(parents=True, exist_ok=True)
    (results / "qc" / "qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
