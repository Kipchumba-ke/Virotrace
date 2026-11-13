"""Gap-mask the alignment for tree inference (never for DRM calling).

A column is masked when more than ``align.mask_gap_frac`` of the QC-passed sequences have a
gap or N there. Masked columns become ``N`` rather than being deleted, so alignment columns
keep their HXB2 coordinates and the tree alignment can be compared with the DRM one.

Output: ``interim/<scope>/aln.pol.masked.fasta`` (QC-passed sequences only) and
``results/<scope>/qc/mask_summary.json``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import read_tsv
from ..settings import Settings

_MISSING = np.zeros(256, dtype=bool)
for _c in "-Nn":
    _MISSING[ord(_c)] = True


def masked_columns(sequences: list[str], max_missing_frac: float) -> np.ndarray:
    """Boolean mask over alignment columns: True where the column is too gappy to use."""
    if not sequences:
        raise ContractError("nothing to mask: no sequences")
    width = {len(s) for s in sequences}
    if len(width) != 1:
        raise ContractError(f"alignment rows have different lengths: {sorted(width)}")
    matrix = np.stack([np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in sequences])
    result: np.ndarray = _MISSING[matrix].mean(axis=0) > max_missing_frac
    return result


def _ranges(keep: np.ndarray) -> list[list[int]]:
    """Inclusive 1-based [start, end] runs of kept (unmasked) columns."""
    out: list[list[int]] = []
    start = None
    for i, k in enumerate([*keep.tolist(), False]):
        if k and start is None:
            start = i
        elif not k and start is not None:
            out.append([start + 1, i])
            start = None
    return out


def run_mask(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    interim = root / settings.paths["interim"] / settings.scope
    results = root / settings.paths["results"] / settings.scope
    aln_path, passed_path = interim / "aln.pol.unmasked.fasta", interim / "qc.passed.tsv"
    for path, stage in ((aln_path, "align"), (passed_path, "qc")):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    passed = {r["strain"] for r in read_tsv(passed_path)}
    rows = [(s, q) for s, q in read_fasta(aln_path) if s in passed]
    if len(rows) != len(passed):
        raise ContractError(f"mask lost records: {len(passed)} passed QC, {len(rows)} aligned")

    gap_frac = float((settings.model_extra or {}).get("align", {}).get("mask_gap_frac", 0.90))
    mask = masked_columns([q for _, q in rows], gap_frac)
    out = np.array(list("N"), dtype="U1")
    masked_rows = []
    for strain, seq in rows:
        chars = np.array(list(seq), dtype="U1")
        chars[mask] = out[0]
        masked_rows.append((strain, "".join(chars.tolist())))
    write_fasta(interim / "aln.pol.masked.fasta", masked_rows, wrap=0)

    summary = {
        "scope": settings.scope,
        "n_sequences": len(rows),
        "n_columns": int(mask.size),
        "n_masked": int(mask.sum()),
        "mask_gap_frac": gap_frac,
        "kept_ranges": _ranges(~mask),
    }
    (results / "qc").mkdir(parents=True, exist_ok=True)
    (results / "qc" / "mask_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
