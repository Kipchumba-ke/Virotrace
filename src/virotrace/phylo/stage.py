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
import re
import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.fasta import read_fasta, write_fasta
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .newick import Tree, midpoint_edge, parse_newick, to_newick
from .temporal import best_root_to_tip, decimal_year, evaluate_gate, permutation_p

Runner = Callable[[Sequence[str]], None]
MIN_SUBTYPE_TIPS = 30


def _excerpt(text: str, head: int = 1500, tail: int = 800) -> str:
    """Start and end of a tool's stderr: the cause is usually at the start, the summary at the end."""
    text = text.strip()
    if len(text) <= head + tail:
        return text
    return f"{text[:head]}\n[...]\n{text[-tail:]}"


def subprocess_runner(command: Sequence[str]) -> None:
    try:
        process = subprocess.run(list(command), capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise ContractError(f"`{command[0]}` is not installed or not on PATH") from exc
    if process.returncode != 0:
        raise ContractError(
            f"`{command[0]}` failed (exit {process.returncode}): {_excerpt(process.stderr)}"
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


def _gate(settings: Settings, tree: Tree, rows: list[dict[str, str]]) -> dict[str, Any]:
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


MIN_PLAUSIBLE_ROOT_YEAR = 1900  # published group M origin estimates are the 1900s to 1930s


def _check_rate_is_plausible(node_data: Path, prior: Sequence[float]) -> None:
    """Reject a time tree whose rate or root date cannot be right.

    A rate several-fold below any HIV-1 pol estimate puts the root centuries back (an A1 tree
    came out at the year 1220), and a rate just inside the prior can still do it (A1 at the year
    1712, G at 1859), so the root date is checked as well: nothing here can be older than group M.
    """
    if not node_data.exists() or len(prior) < 2:
        return
    data = json.loads(node_data.read_text())
    rate = (data.get("clock") or {}).get("rate")
    if rate is not None:
        low, high = min(prior) / 5, max(prior) * 5
        if not low <= rate <= high:
            raise ContractError(
                f"estimated clock rate {rate:.3g} substitutions/site/year is outside "
                f"{low:.3g} to {high:.3g} (the configured prior range, widened five-fold)"
            )
    dates = [n["numdate"] for n in data.get("nodes", {}).values() if "numdate" in n]
    if dates and min(dates) < MIN_PLAUSIBLE_ROOT_YEAR:
        raise ContractError(
            f"the estimated root date {min(dates):.0f} is before {MIN_PLAUSIBLE_ROOT_YEAR}, "
            "older than HIV-1 group M can be"
        )


def run_refine(
    settings: Settings,
    *,
    repo_root: Path | str = ".",
    runner: Runner = subprocess_runner,
    subtype: str | None = None,
    out_dir: Path | None = None,
    gate: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """``out_dir`` writes the refine outputs elsewhere (a sensitivity run); ``gate`` reuses an
    existing gate result instead of recomputing it, permutation test included."""
    p = _paths(settings, Path(repo_root), subtype)
    for key, stage in (("tree", "phylo tree"), ("meta", "subsample"), ("aln", "subsample")):
        if not p[key].exists():
            raise ContractError(f"{p[key]} not found; run `virotrace {stage}` first")
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        p = {**p, "phylo": out_dir}
    rows = read_tsv(p["meta"])
    tree = parse_newick(p["tree"].read_text())
    gate_cfg = settings.phylo.temporal_gate
    if gate is not None:
        gate = dict(gate)
    elif gate_cfg.enabled:
        gate = _gate(settings, tree, rows)
    else:
        gate = {"pooled": {"passed": True}}
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
    rooted = p["phylo"] / "tree.midpoint.nwk"

    def midpoint_rooted(cmd: list[str]) -> list[str]:
        """The same command on the midpoint-rooted tree, kept as rooted."""
        rooted.write_text(to_newick(tree, midpoint_edge(tree)) + "\n")
        out = list(cmd)
        out[out.index("--tree") + 1] = str(rooted)
        return [*out, "--keep-root"]

    if use_time:
        command += [
            "--timetree", "--root", "best", "--coalescent", "opt",
            "--date-inference", "marginal", "--clock-filter-iqd", str(cfg.clock_filter_iqd),
        ]  # fmt: skip
        if cfg.clock_rate is not None:
            command += ["--clock-rate", str(cfg.clock_rate)]
        if cfg.date_confidence:
            command.append("--date-confidence")
        check = cfg.clock_rate is None

        def attempt(cmd: list[str]) -> None:
            node_data.unlink(missing_ok=True)  # never judge a rate left by an earlier run
            runner(cmd)
            if check:
                _check_rate_is_plausible(node_data, cfg.clock_rate_prior)

        try:
            attempt(command)
        except ContractError as first:
            # TreeTime can fail while searching for the best root (e.g. "Rerooting failed")
            # even where the gate's own regression is fine. Try the same time tree on the
            # midpoint root; failing that, give a divergence tree. Either is said in the gate.
            gate["refine_error"] = str(first)[:2500]
            no_best_root = [
                c for i, c in enumerate(command) if not (c == "best" and command[i - 1] == "--root")
            ]
            no_best_root.remove("--root")
            retry = midpoint_rooted(no_best_root)
            try:
                attempt(retry)
                command = retry
                gate["root"] = "midpoint root; augur's best-root search failed"
                gate["note"] = "time tree rooted at the midpoint: augur's best-root search failed"
            except ContractError as second:
                gate["refine_error_midpoint"] = str(second)[:2500]
                # TreeTime's own rate estimate (after its clock filter) can be negative even
                # where the gate's regression is positive; it then asks for the rate and its
                # standard deviation. Give it the gate's regression slope, with a wide spread.
                slope = (gate.get("pooled") or {}).get("slope")
                if isinstance(slope, (int, float)) and slope > 0:
                    fixed = [
                        *retry,
                        "--clock-rate",
                        f"{slope:.6g}",
                        "--clock-std-dev",
                        f"{slope / 2:.6g}",
                    ]
                    if "--clock-rate" in retry:
                        i = retry.index("--clock-rate")
                        fixed = [c for j, c in enumerate(retry) if j not in (i, i + 1)] + fixed[
                            len(retry) :
                        ]
                    try:
                        runner(fixed)
                        command = fixed
                        gate["root"] = (
                            "midpoint root; augur's own rate estimate failed or was implausible"
                        )
                        gate["fixed_clock_rate"] = slope
                        gate["note"] = (
                            "time tree at the midpoint root with the clock rate fixed at the "
                            "regression slope (std dev half of it): augur's own estimate failed"
                        )
                    except ContractError as third:
                        gate["refine_error_fixed_rate"] = str(third)[:2500]
                if "fixed_clock_rate" not in gate:
                    gate["time_tree"] = False
                    gate["note"] = "augur could not date this tree; showing a divergence tree"
                    command = midpoint_rooted([c for c in command[: command.index("--timetree")]])
                    runner(command)
    else:
        # no clock to root on: midpoint root, written out so augur keeps it as given
        command = midpoint_rooted(command)
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
        try:
            run_tree(settings, repo_root=repo_root, runner=runner, subtype=subtype)
            result = run_refine(settings, repo_root=repo_root, runner=runner, subtype=subtype)
        except ContractError as exc:
            raise ContractError(f"subtype {subtype}: {exc}") from exc
        gate = result["temporal_gate"]
        out[subtype] = {
            "time_tree": gate["time_tree"],
            "r_squared": gate["pooled"].get("r_squared"),
            "permutation_p": gate["pooled"].get("permutation_p"),
            "slope": gate["pooled"].get("slope"),
        }
    return {"scope": settings.scope, "subtype_trees": out}


def _clock_summary(node_data: Path) -> dict[str, float | None]:
    data = json.loads(node_data.read_text())
    dates = [n["numdate"] for n in data.get("nodes", {}).values() if "numdate" in n]
    return {
        "clock_rate": (data.get("clock") or {}).get("rate"),
        "root_date": round(min(dates), 2) if dates else None,
        "latest_date": round(max(dates), 2) if dates else None,
    }


def run_clock_sensitivity(
    settings: Settings,
    rates: Sequence[float] | None = None,
    *,
    repo_root: Path | str = ".",
    runner: Runner = subprocess_runner,
    subtypes: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Refit each dated tree with the clock rate fixed, to see how much its dates depend on it.

    Estimated rates in this corpus are low (A1 about 0.0005 and D about 0.0007 substitutions per
    site per year, below the configured prior of 0.0008 to 0.0030), so the root date is the
    quantity that moves. The default rates are the two ends of ``refine.clock_rate_prior``. A tree that
    failed its temporal gate is refused: fixing a rate on it would manufacture a timeline. The
    topology is reused; only the dating is redone, into ``clock_<rate>/`` beside the main run.
    """
    cfg = settings.phylo.refine
    rates = list(rates) if rates else list(cfg.clock_rate_prior)
    targets: list[str | None] = list(subtypes) if subtypes else list(cfg.subtype_trees) or [None]
    out: dict[str, Any] = {}
    for subtype in targets:
        p = _paths(settings, Path(repo_root), subtype)
        gate_path, base_data = p["phylo"] / "temporal_gate.json", p["phylo"] / "branch_lengths.json"
        label = subtype or "pooled"
        step = "phylo refine" if subtype is None else "phylo subtrees"
        if not gate_path.exists():
            raise ContractError(f"{p['phylo']} has no refined tree; run `virotrace {step}` first")
        gate = json.loads(gate_path.read_text())
        if not gate.get("time_tree"):
            raise ContractError(
                f"the {label} tree failed its temporal gate, so it has no dates to test; "
                "a fixed rate would not repair that"
            )
        if not base_data.exists():
            raise ContractError(f"{base_data} not found; run `virotrace {step}` first")
        base = _clock_summary(base_data)
        if gate.get("fixed_clock_rate"):
            base["rate_fixed"] = True  # the pipeline supplied this rate; it was not estimated
        runs = {}
        for rate in rates:
            fixed = settings.model_copy(
                update={
                    "phylo": settings.phylo.model_copy(
                        update={"refine": cfg.model_copy(update={"clock_rate": rate})}
                    )
                }
            )
            folder = p["phylo"] / f"clock_{rate:g}"
            run_refine(
                fixed, repo_root=repo_root, runner=runner, subtype=subtype,
                out_dir=folder, gate=gate,
            )  # fmt: skip
            summary = _clock_summary(folder / "branch_lengths.json")
            shift = (
                round(summary["root_date"] - base["root_date"], 2)
                if summary["root_date"] is not None and base["root_date"] is not None
                else None
            )
            runs[f"{rate:g}"] = {**summary, "root_shift_years": shift}
        out[label] = {"estimated": base, "fixed": runs}
        result = {"scope": settings.scope, "tree": label, **out[label]}
        (p["phylo"] / "clock_sensitivity.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
    return {"scope": settings.scope, "clock_sensitivity": out}


def run_subtype_tree(
    settings: Settings, *, repo_root: Path | str = ".", runner: Runner = subprocess_runner
) -> dict[str, Any]:
    """Second subtype call from a tree of the pure-subtype panel plus the tree set."""
    from ..curate.subtype import agreement, panel_label, parent_group, pure_rows
    from .subtype_tree import concordance, tree_subtype_calls

    root = Path(repo_root)
    p = _paths(settings, root)
    panel_path = root / settings.paths["interim"] / settings.scope / "panel.pol.fasta"
    for path, stage in ((panel_path, "subtype"), (p["aln"], "subsample"), (p["meta"], "subsample")):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")

    query = list(read_fasta(p["aln"]))
    unusable = [all(c in "N-" for c in col) for col in zip(*(q for _, q in query), strict=True)]
    panel_rows = pure_rows(list(read_fasta(panel_path)))
    if len(panel_rows) < 4:
        raise ContractError("the pure-subtype panel has too few references to build a tree")
    width = len(query[0][1])

    def safe(n: str) -> str:
        return "ref_" + re.sub(r"[^A-Za-z0-9_.\-]", "_", n)

    labelled: dict[str, str] = {}
    rows: list[tuple[str, str]] = []
    for name, seq in panel_rows:
        if len(seq) != width:
            raise ContractError(f"panel width {len(seq)} != tree alignment width {width}")
        masked = "".join("N" if u else c for c, u in zip(seq, unusable, strict=True))
        labelled[safe(name)] = panel_label(name)
        rows.append((safe(name), masked))
    out_dir = p["phylo"] / "subtype_tree"
    out_dir.mkdir(parents=True, exist_ok=True)
    aln = out_dir / "aln.fasta"
    write_fasta(aln, rows + query, wrap=0)

    prefix = out_dir / "iqtree"
    for stale in out_dir.glob("iqtree.*"):
        stale.unlink()
    runner(iqtree_command(settings, aln, prefix))
    treefile = prefix.with_suffix(".treefile")
    if not treefile.exists():
        raise ContractError(f"{treefile} was not written; see {prefix}.log")
    tree = parse_newick(treefile.read_text())

    groups = {lab: parent_group(lab) for lab in set(labelled.values())}
    cfg = settings.subtype
    calls = tree_subtype_calls(tree, labelled, groups, min_support=cfg.tree_min_support)
    meta = {r["strain"]: r for r in read_tsv(p["meta"])}
    table = []
    for strain, c in calls.items():
        m = meta.get(strain, {})
        table.append(
            {
                "strain": strain, "subtype_tree": c["call"], "tree_level": c["level"],
                "tree_support": "" if c["support"] is None else c["support"],
                "subtype_called": m.get("subtype_called", ""), "crf_label": m.get("crf_label", ""),
                "subtype_reported": m.get("subtype_reported", ""),
            }
        )  # fmt: skip
    summary, confusion = concordance(table, agreement)
    write_tsv(out_dir / "subtype_calls_tree.tsv", table)
    write_tsv(
        out_dir / "subtype_disagreements.tsv",
        [r for r in table if r["concordance"] == "disagree"],
        list(table[0]) if table else None,
    )
    result = {
        "scope": settings.scope,
        "n_queries": len(table),
        "n_references": len(rows),
        "concordance": summary,
        "confusion_called_vs_tree": confusion,
        "min_support": cfg.tree_min_support,
    }
    (out_dir / "subtype_tree_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result
