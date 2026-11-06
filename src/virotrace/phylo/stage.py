"""Tree inference and the temporal-signal gate.

``phylo tree``    IQ-TREE on ``data/processed/<scope>/aln.tree.fasta``.
``phylo refine``  root-to-tip regression with the best root (the gate), then ``augur refine``:
                  a time tree if the gate passes, a divergence tree with the reason recorded
                  if it does not.

External programs are injected (``Runner``), so everything around them is testable without
IQ-TREE or augur installed.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .newick import midpoint_edge, parse_newick, to_newick
from .temporal import best_root_to_tip, decimal_year, evaluate_gate, permutation_p

Runner = Callable[[Sequence[str]], None]
MIN_SUBTYPE_TIPS = 30


def subprocess_runner(command: Sequence[str]) -> None:
    try:
        process = subprocess.run(list(command), capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise ContractError(f"`{command[0]}` is not installed or not on PATH") from exc
    if process.returncode != 0:
        raise ContractError(
            f"`{command[0]}` failed (exit {process.returncode}): {process.stderr.strip()[-800:]}"
        )


def iqtree_command(settings: Settings, aln: Path, prefix: Path) -> list[str]:
    cfg = settings.phylo.tree
    command = [
        cfg.tool,
        "-s", str(aln),
        "-m", "MFP" if cfg.model_finder else cfg.model,
        "-T", str(settings.threads),
        "--seed", str(settings.stage_seed("tree")),
        "--prefix", str(prefix),
        "-redo",
    ]  # fmt: skip
    if cfg.bootstrap == "ufboot":
        command += ["-B", str(cfg.replicates)]
    if cfg.fast_mode:
        command.append("--fast")
    return command


def _paths(settings: Settings, root: Path, subtype: str | None = None) -> dict[str, Path]:
    processed = root / settings.paths["processed"] / settings.scope
    phylo = root / settings.paths["results"] / settings.scope / "phylo"
    if subtype is None:
        return {
            "aln": processed / "aln.tree.fasta",
            "meta": processed / "metadata.final.tsv",
            "phylo": phylo,
            "tree": phylo / "tree.nwk",
        }
    return {
        "aln": processed / f"aln.tree.{subtype}.fasta",
        "meta": processed / f"metadata.final.{subtype}.tsv",
        "phylo": phylo / subtype,
        "tree": phylo / subtype / "tree.nwk",
    }


def _write_subtype_inputs(settings: Settings, root: Path, subtype: str) -> int:
    """Subset the tree alignment and metadata to one called subtype. Returns the tip count."""
    base, sub = _paths(settings, root), _paths(settings, root, subtype)
    for key in ("aln", "meta"):
        if not base[key].exists():
            raise ContractError(f"{base[key]} not found; run `virotrace subsample` first")
    rows = [r for r in read_tsv(base["meta"]) if r.get("subtype_called") == subtype]
    if len(rows) < 3:
        raise ContractError(f"only {len(rows)} tips called {subtype}; cannot build a tree")
    keep = {r["strain"] for r in rows}
    write_fasta(sub["aln"], [(n, q) for n, q in read_fasta(base["aln"]) if n in keep], wrap=0)
    write_tsv(sub["meta"], rows)
    return len(rows)


def run_tree(
    settings: Settings,
    *,
    repo_root: Path | str = ".",
    runner: Runner = subprocess_runner,
    subtype: str | None = None,
) -> dict[str, Any]:
    if subtype is not None:
        _write_subtype_inputs(settings, Path(repo_root), subtype)
    p = _paths(settings, Path(repo_root), subtype)
    if not p["aln"].exists():
        raise ContractError(f"{p['aln']} not found; run `virotrace subsample` first")
    p["phylo"].mkdir(parents=True, exist_ok=True)
    prefix = p["phylo"] / "iqtree"
    command = iqtree_command(settings, p["aln"], prefix)
    for stale in prefix.parent.glob(f"{prefix.name}.*"):  # never read a previous run's tree
        stale.unlink()
    runner(command)
    treefile = prefix.with_suffix(".treefile")
    if not treefile.exists():
        raise ContractError(f"IQ-TREE produced no tree at {treefile}")
    shutil.copyfile(treefile, p["tree"])
    parse_newick(p["tree"].read_text())  # fail here, not three stages later
    return {
        "scope": settings.scope, "subtype": subtype, "command": command, "tree": str(p["tree"])
    }  # fmt: skip


def augur_date(date: str, precision: str) -> str:
    """Augur's date format: ambiguous parts as XX."""
    parts = date.split("-") if date else []
    if precision == "day" and len(parts) == 3:
        return date
    if precision == "month" and len(parts) >= 2:
        return f"{parts[0]}-{parts[1]}-XX"
    if precision == "year" and parts:
        return f"{parts[0]}-XX-XX"
    return "?"


def _gate(settings: Settings, tree, rows: list[dict[str, str]]) -> dict[str, Any]:
    by_strain = {r["strain"]: r for r in rows}
    missing = [n for n in (tree.name[t] for t in tree.tips) if n not in by_strain]
    if missing:
        raise ContractError(
            f"{len(missing)} tree tips are not in metadata.final.tsv, e.g. {missing[:3]}"
        )
    tip_rows = [by_strain[tree.name[t]] for t in tree.tips]
    dates = [decimal_year(r["date"]) for r in tip_rows]
    usable = [d is not None for d in dates]
    gate_cfg = settings.phylo.temporal_gate
    y = [d if d is not None else 0.0 for d in dates]

    def run(mask: list[bool]) -> dict[str, Any]:
        used, d = best_root_to_tip(
            tree, y, subset=mask, require_positive_slope=gate_cfg.require_positive_slope
        )
        result = evaluate_gate(
            used, d, min_r2=gate_cfg.min_temporal_r2,
            require_positive_slope=gate_cfg.require_positive_slope,
        )  # fmt: skip
        out = {
            "n": result.n, "slope": result.slope, "r_squared": round(result.r_squared, 4),
            "passed": result.passed, "reason": result.reason,
        }  # fmt: skip
        if gate_cfg.date_randomisation:
            p = permutation_p(
                tree, y, result.r_squared, subset=mask, n=gate_cfg.n_permutations,
                seed=settings.stage_seed("temporal_gate"),
                require_positive_slope=gate_cfg.require_positive_slope,
            )  # fmt: skip
            out["permutation_p"] = round(p, 4)
            if result.passed and p > gate_cfg.max_permutation_p:
                out["passed"] = False
                out["reason"] = (
                    f"R^2 {result.r_squared:.3f} is not distinguishable from shuffled dates "
                    f"(permutation p = {p:.3f}, limit {gate_cfg.max_permutation_p})"
                )
        return out

    out = {
        "n_tips": len(tip_rows),
        "n_undated": usable.count(False),
        "undated_examples": [
            f"{r['strain']}: date={r['date']!r} precision={r['date_precision']!r}"
            for r, ok in zip(tip_rows, usable, strict=True)
            if not ok
        ][:5],
        "min_r2": gate_cfg.min_temporal_r2,
        "pooled": run(usable),
    }
    # a pooled tree is dominated by between-subtype divergence, so also ask within each subtype
    by_subtype = {}
    subtypes = [r.get("subtype_called", "") for r in tip_rows]
    for label in sorted(set(subtypes) - {"", "unassigned"}):
        mask = [u and s == label for u, s in zip(usable, subtypes, strict=True)]
        if sum(mask) >= MIN_SUBTYPE_TIPS and len({y[i] for i, m in enumerate(mask) if m}) > 1:
            by_subtype[label] = run(mask)
    out["by_subtype"] = by_subtype
    return out


def run_refine(
    settings: Settings,
    *,
    repo_root: Path | str = ".",
    runner: Runner = subprocess_runner,
    subtype: str | None = None,
) -> dict[str, Any]:
    p = _paths(settings, Path(repo_root), subtype)
    for key, stage in (("tree", "phylo tree"), ("meta", "subsample"), ("aln", "subsample")):
        if not p[key].exists():
            raise ContractError(f"{p[key]} not found; run `virotrace {stage}` first")
    rows = read_tsv(p["meta"])
    tree = parse_newick(p["tree"].read_text())
    gate_cfg = settings.phylo.temporal_gate
    gate = _gate(settings, tree, rows) if gate_cfg.enabled else {"pooled": {"passed": True}}
    use_time = settings.phylo.refine.enabled and (not gate_cfg.enabled or gate["pooled"]["passed"])
    gate["time_tree"] = use_time
    if gate_cfg.enabled and not gate["pooled"]["passed"]:
        gate["note"] = gate["pooled"]["reason"]

    write_tsv(
        p["phylo"] / "metadata.augur.tsv",
        [{"strain": r["strain"], "date": augur_date(r["date"], r["date_precision"])} for r in rows],
        ["strain", "date"],
    )
    refined, node_data = p["phylo"] / "tree.refined.nwk", p["phylo"] / "branch_lengths.json"
    command = [
        "augur", "refine",
        "--tree", str(p["tree"]), "--alignment", str(p["aln"]),
        "--metadata", str(p["phylo"] / "metadata.augur.tsv"), "--metadata-id-columns", "strain",
        "--output-tree", str(refined), "--output-node-data", str(node_data),
        "--seed", str(settings.stage_seed("refine")),
    ]  # fmt: skip
    cfg = settings.phylo.refine
    if use_time:
        command += [
            "--timetree", "--root", "best", "--coalescent", "opt",
            "--date-inference", "marginal", "--clock-filter-iqd", str(cfg.clock_filter_iqd),
        ]  # fmt: skip
        if cfg.clock_rate is not None:
            command += ["--clock-rate", str(cfg.clock_rate)]
        if cfg.date_confidence:
            command.append("--date-confidence")
    else:
        # no clock to root on: midpoint root, written out so augur keeps it as given
        rooted = p["phylo"] / "tree.midpoint.nwk"
        rooted.write_text(to_newick(tree, midpoint_edge(tree)) + "\n")
        command[command.index("--tree") + 1] = str(rooted)
        command.append("--keep-root")
    runner(command)
    (p["phylo"] / "temporal_gate.json").write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n"
    )
    return {"scope": settings.scope, "subtype": subtype, "temporal_gate": gate, "command": command}


def run_subtype_trees(
    settings: Settings, *, repo_root: Path | str = ".", runner: Runner = subprocess_runner
) -> dict[str, Any]:
    """One tree per configured subtype, each with its own gate (and time tree if it passes)."""
    out: dict[str, Any] = {}
    for subtype in settings.phylo.refine.subtype_trees:
        run_tree(settings, repo_root=repo_root, runner=runner, subtype=subtype)
        result = run_refine(settings, repo_root=repo_root, runner=runner, subtype=subtype)
        gate = result["temporal_gate"]
        out[subtype] = {
            "time_tree": gate["time_tree"],
            "r_squared": gate["pooled"].get("r_squared"),
            "permutation_p": gate["pooled"].get("permutation_p"),
            "slope": gate["pooled"].get("slope"),
        }
    return {"scope": settings.scope, "subtype_trees": out}
