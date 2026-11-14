"""Ancestral trait reconstruction (``augur traits``) over country and region.

Runs on the pooled refined tree. Subnational (admin1) reconstruction is added only if enough
sequences have a resolved admin1 (``phylo.traits.admin1``); otherwise the geographic layer
stops at country, and the summary records the resolved fraction so the report can say so.

Discrete-trait reconstruction on a convenience sample partly reflects where sequences were
deposited, not where the virus moved. The result is descriptive and the limitation is stated
in the report.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .stage import Runner, subprocess_runner

TRUE = {"true", "1", "yes"}


def resolved_fraction(rows: list[dict[str, str]]) -> float:
    """Share of rows with a resolved admin1. Zero when there are no rows."""
    if not rows:
        return 0.0
    return sum(r.get("admin1_resolved", "").strip().lower() in TRUE for r in rows) / len(rows)


def trait_columns(
    settings: Settings, rows: list[dict[str, str]]
) -> tuple[list[str], dict[str, Any]]:
    cfg = settings.phylo.traits
    fraction = resolved_fraction(rows)
    mode = cfg.admin1.enabled
    use = mode is True or (mode == "auto" and fraction >= cfg.admin1.min_resolved_frac)
    columns = list(cfg.columns)
    if use and "admin1" not in columns:
        columns.append("admin1")
    return columns, {
        "mode": mode,
        "resolved_fraction": round(fraction, 4),
        "min_resolved_frac": cfg.admin1.min_resolved_frac,
        "reconstructed": use,
    }


def run_traits(
    settings: Settings,
    *,
    repo_root: Path | str = ".",
    runner: Runner = subprocess_runner,
) -> dict[str, Any]:
    root = Path(repo_root)
    phylo = root / settings.paths["results"] / settings.scope / "phylo"
    tree = phylo / "tree.refined.nwk"
    meta = root / settings.paths["processed"] / settings.scope / "metadata.final.tsv"
    for path, stage in ((tree, "phylo refine"), (meta, "subsample")):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    rows = read_tsv(meta)
    columns, admin1 = trait_columns(settings, rows)
    absent = [c for c in columns if rows and c not in rows[0]]
    if absent:
        raise ContractError(f"metadata.final.tsv has no column(s) {absent} to reconstruct")

    phylo.mkdir(parents=True, exist_ok=True)
    table = phylo / "metadata.traits.tsv"
    write_tsv(
        table,
        [
            {"strain": r["strain"], **{c: (r.get(c) or "").strip() or "?" for c in columns}}
            for r in rows
        ],
        ["strain", *columns],
    )
    output = phylo / "traits.json"
    output.unlink(missing_ok=True)
    command = [
        "augur", "traits",
        "--tree", str(tree), "--metadata", str(table), "--metadata-id-columns", "strain",
        "--columns", *columns, "--confidence", "--output-node-data", str(output),
    ]  # fmt: skip
    runner(command)
    if not output.exists():
        raise ContractError(f"augur traits produced no node data at {output}")
    json.loads(output.read_text())  # fail here, not at export

    summary = {
        "scope": settings.scope,
        "columns": columns,
        "n_tips": len(rows),
        "n_missing": {c: sum(not (r.get(c) or "").strip() for r in rows) for c in columns},
        "admin1": admin1,
        "traits": str(output),
    }
    (phylo / "traits_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
