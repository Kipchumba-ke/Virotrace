"""The align stage: dedup'd sequences -> pol alignment in HXB2 coordinates.

Outputs (per scope):
    interim/<scope>/aln.pol.unmasked.fasta   one row per sequence, one column per HXB2
                                             pol position (3012); never masked
    interim/<scope>/insertions.tsv           insertions relative to HXB2 inside pol
    results/<scope>/qc/align_summary.json    counts, including record conservation

The alignment is deliberately *not* masked here. Masking is an acceptable
tree-building convenience and an unacceptable resistance-calling shortcut, so the
masked tree alignment is produced later from this one.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import write_tsv
from ..settings import Settings
from .mafft import Runner, align_chunk, chunked, mafft_runner
from .project import project_to_reference

POL_GENE = "pol"


def run_align(
    settings: Settings,
    *,
    repo_root: Path | str = ".",
    runner: Runner = mafft_runner,
) -> dict[str, Any]:
    from ..drm.reference import load_reference_sequence

    root = Path(repo_root)
    scope = settings.scope
    pol_start, pol_end = settings.reference.genes.pol
    chunk_size = int((settings.model_extra or {}).get("align", {}).get("chunk_size", 500))

    interim = root / settings.paths["interim"] / scope
    source = interim / "seqs.dedup.fasta"
    if not source.exists():
        raise ContractError(f"{source} not found; run `virotrace normalise` first")

    record = root / settings.paths["reference"] / settings.reference.record
    reference_pol = load_reference_sequence(record)[pol_start - 1 : pol_end]
    if len(reference_pol) != pol_end - pol_start + 1:
        raise ContractError(
            f"{record} is shorter than the configured pol span {pol_start}..{pol_end}"
        )
    reference = (f"{settings.reference.accession}_pol", reference_pol)

    projected_rows: list[tuple[str, str]] = []
    insertion_rows: list[dict[str, Any]] = []
    n_in = n_chunks = 0

    for chunk in chunked(read_fasta(source), chunk_size):
        n_chunks += 1
        n_in += len(chunk)
        aligned_ref, rows = align_chunk(reference, chunk, threads=settings.threads, runner=runner)
        for strain, _ in chunk:  # input order, so the output is deterministic
            projection = project_to_reference(
                aligned_ref, rows[strain], start=1, end=len(reference_pol)
            )
            projected_rows.append((strain, projection.projected))
            insertion_rows.extend(
                {"strain": strain, "after_ref_pos": ins.after_ref_pos, "length": ins.length}
                for ins in projection.insertions
            )

    if n_in == 0:
        raise ContractError(f"{source} contains no sequences")
    if len(projected_rows) != n_in:
        raise ContractError(f"alignment lost records: {n_in} in, {len(projected_rows)} out")

    write_fasta(interim / "aln.pol.unmasked.fasta", projected_rows, wrap=0)
    write_tsv(interim / "insertions.tsv", insertion_rows, ["strain", "after_ref_pos", "length"])

    summary = {
        "scope": scope,
        "reference": reference[0],
        "pol_span": [pol_start, pol_end],
        "chunk_size": chunk_size,
        "n_chunks": n_chunks,
        "n_sequences": n_in,
        "n_with_insertions": len({row["strain"] for row in insertion_rows}),
    }
    qc_dir = root / settings.paths["results"] / scope / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    (qc_dir / "align_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
