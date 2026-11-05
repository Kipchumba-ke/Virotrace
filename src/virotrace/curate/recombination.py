"""Recombination screening by sliding-window subtype scanning over pure-subtype parent groups.

For each QC-passed sequence, windows of ``window`` alignment columns (step ``step``) are
placed against the subtype panel. A window is *decisive* when its nearest label beats the
next label by at least ``min_margin`` (p-distance). Consecutive decisive windows with the
same label form a segment; a segment must span ``min_segment`` columns to count.

* ``suspected``  -- two or more supported segments with different labels
* ``confirmed``  -- as above, and each segment, re-placed on its own columns with site
                    bootstrap, keeps its label at ``recombination.confirm_support``
* ``none``       -- one supported label across the sequence
* ``not_testable`` -- the sequence spans fewer than 2 x ``min_segment`` columns, so two
                      supported segments are impossible

This is a screen with real limits, stated in the report: pol is short (a PR+RT genotype is
~1.3 kb), so recent or intra-subtype recombination and breakpoints near the ends are
invisible, and a mosaic of lineages absent from the panel looks like noise, not like a
recombinant. It is not 3SEQ or RDP and does not claim their statistics.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from ..errors import ContractError
from ..io.fasta import read_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .subtype import Panel, build_panel, encode, panel_label, parent_group, place, pure_rows

FLAG_COLUMNS = [
    "strain",
    "recomb_status",
    "n_windows",
    "n_decisive",
    "segments",
    "exclude_from_tree",
]


def scan_windows(
    query: np.ndarray, panel: Panel, *, window: int, step: int, min_sites: int, min_margin: float
) -> list[tuple[int, int, str | None]]:
    """``(start, end, label)`` per testable window; label is ``None`` when not decisive."""
    out: list[tuple[int, int, str | None]] = []
    rng = np.random.default_rng(0)  # no bootstrap here, so the stream is never consumed
    for start in range(0, panel.length - window + 1, step):
        columns = np.arange(start, start + window)
        placed = place(query, panel, n_bootstrap=0, min_sites=min_sites, rng=rng, columns=columns)
        if placed is None:
            continue
        decisive = (
            placed.runner_up == "" or placed.runner_up_distance - placed.distance >= min_margin
        )
        out.append((start, start + window, placed.label if decisive else None))
    return out


def segments_from_windows(
    windows: list[tuple[int, int, str | None]], *, min_segment: int
) -> list[tuple[str, int, int]]:
    """Merge consecutive decisive windows of one label; keep those spanning ``min_segment``."""
    runs: list[list[Any]] = []
    for start, end, label in windows:
        if label is None:
            continue
        if runs and runs[-1][0] == label:
            runs[-1][2] = end
        else:
            runs.append([label, start, end])
    return [(label, s, e) for label, s, e in runs if e - s >= min_segment]


def classify(
    query: np.ndarray, panel: Panel, cfg: Any, *, rng: np.random.Generator
) -> tuple[str, list[tuple[str, int, int]], int, int]:
    covered = np.flatnonzero(query < 4)
    if covered.size == 0 or covered[-1] - covered[0] + 1 < 2 * cfg.min_segment:
        return "not_testable", [], 0, 0  # cannot hold two supported segments
    windows = scan_windows(
        query,
        panel,
        window=cfg.window,
        step=cfg.step,
        min_sites=cfg.min_window_sites,
        min_margin=cfg.min_margin,
    )
    if len(windows) < 2:
        return "not_testable", [], len(windows), 0
    decisive = sum(1 for _, _, label in windows if label is not None)
    segments = segments_from_windows(windows, min_segment=cfg.min_segment)
    # collapse adjacent segments that kept the same label after short interruptions
    merged: list[tuple[str, int, int]] = []
    for label, s, e in segments:
        if merged and merged[-1][0] == label:
            merged[-1] = (label, merged[-1][1], e)
        else:
            merged.append((label, s, e))
    if len({label for label, _, _ in merged}) < 2:
        return "none", merged, len(windows), decisive

    for label, s, e in merged:
        placed = place(
            query,
            panel,
            n_bootstrap=cfg.confirm_bootstrap,
            min_sites=cfg.min_window_sites,
            rng=rng,
            columns=np.arange(s, e),
        )
        if placed is None or placed.label != label or placed.support < cfg.confirm_support:
            return "suspected", merged, len(windows), decisive
    return "confirmed", merged, len(windows), decisive


def excluded_from_tree(status: str, policy: str) -> bool:
    return (policy == "exclude_confirmed" and status == "confirmed") or (
        policy == "exclude_suspected" and status in ("suspected", "confirmed")
    )


def run_recombination(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    interim = root / settings.paths["interim"] / scope
    results = root / settings.paths["results"] / scope
    cfg = settings.recombination
    for path, stage in (
        (interim / "aln.pol.unmasked.fasta", "align"),
        (interim / "panel.pol.fasta", "subtype"),
        (interim / "subtype_calls.tsv", "subtype"),
    ):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    if not cfg.enabled:
        status_rows: list[dict[str, Any]] = [
            {
                "strain": r["strain"],
                "recomb_status": "not_run",
                "n_windows": 0,
                "n_decisive": 0,
                "segments": "",
                "exclude_from_tree": False,
            }
            for r in read_tsv(interim / "subtype_calls.tsv")
        ]
        write_tsv(interim / "recombination.tsv", status_rows, FLAG_COLUMNS)
        return {"scope": scope, "enabled": False, "n_sequences": len(status_rows)}

    # Parents are the pure-subtype groups (A1..A8 -> A, F1/F2 -> F). CRFs are not parents:
    # a CRF is itself a mosaic of them, and listing it as a parent makes windows non-decisive.
    panel = build_panel(
        [
            (f"{parent_group(panel_label(name))}.{name}", seq)
            for name, seq in pure_rows(list(read_fasta(interim / "panel.pol.fasta")))
        ]
    )
    called = {r["strain"] for r in read_tsv(interim / "subtype_calls.tsv")}
    rng = np.random.default_rng(settings.stage_seed("recombination"))

    rows: list[dict[str, Any]] = []
    for strain, sequence in read_fasta(interim / "aln.pol.unmasked.fasta"):
        if strain not in called:
            continue
        status, segments, n_windows, n_decisive = classify(encode(sequence), panel, cfg, rng=rng)
        rows.append(
            {
                "strain": strain,
                "recomb_status": status,
                "n_windows": n_windows,
                "n_decisive": n_decisive,
                "segments": ";".join(f"{label}:{s}-{e}" for label, s, e in segments),
                "exclude_from_tree": excluded_from_tree(status, cfg.tree_policy),
            }
        )
    if len(rows) != len(called):
        raise ContractError(f"recombination lost records: {len(called)} called, {len(rows)} out")
    write_tsv(interim / "recombination.tsv", rows, FLAG_COLUMNS)

    summary = {
        "scope": scope,
        "n_sequences": len(rows),
        "status": dict(Counter(r["recomb_status"] for r in rows)),
        # segment order along pol for flagged sequences, e.g. "A>D": count
        "segment_patterns": dict(
            Counter(
                ">".join(seg.split(":")[0] for seg in r["segments"].split(";"))
                for r in rows
                if r["recomb_status"] in ("confirmed", "suspected")
            ).most_common(20)
        ),
        "tree_policy": cfg.tree_policy,
        "drm_policy": cfg.drm_policy,
        "n_excluded_from_tree": sum(1 for r in rows if r["exclude_from_tree"]),
        "window": cfg.window,
        "step": cfg.step,
        "min_segment": cfg.min_segment,
    }
    (results / "qc").mkdir(parents=True, exist_ok=True)
    (results / "qc" / "recombination_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
