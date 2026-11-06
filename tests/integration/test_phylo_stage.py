"""Newick, root-to-tip with the best root, and the tree/refine stages with injected runners."""

from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import pytest

from virotrace.errors import ContractError
from virotrace.io.fasta import write_fasta
from virotrace.io.tables import write_tsv
from virotrace.phylo.newick import midpoint_edge, parse_newick, tip_distances, to_newick
from virotrace.phylo.stage import augur_date, iqtree_command, run_refine, run_tree
from virotrace.phylo.temporal import (
    best_root_to_tip,
    decimal_year,
    evaluate_gate,
    permutation_p,
)
from virotrace.settings import load_settings


def test_parse_newick_ignores_support_labels_and_numbers_tips_in_preorder():
    tree = parse_newick("((A:1,B:2)95.5/100:0.5,C:3,D:4);")
    assert [tree.name[t] for t in tree.tips] == ["A", "B", "C", "D"]
    assert tree.length[tree.tips[1]] == 2.0
    with pytest.raises(ContractError):
        parse_newick("(A:1,B:2")  # no terminator
    with pytest.raises(ContractError):
        parse_newick("(A:1,B:2);")  # fewer than three tips


def _random_newick(rng: random.Random, names: list[str]) -> str:
    if len(names) == 1:
        return f"{names[0]}:{rng.uniform(0.01, 1):.4f}"
    cut = rng.randint(1, len(names) - 1)
    left, right = _random_newick(rng, names[:cut]), _random_newick(rng, names[cut:])
    return f"({left},{right}):{rng.uniform(0.01, 1):.4f}"


def test_tip_distances_match_brute_force_path_lengths():
    rng = random.Random(3)
    names = [f"t{i}" for i in range(9)]
    tree = parse_newick(f"({_random_newick(rng, names[:4])},{_random_newick(rng, names[4:])});")
    dist, lo, hi = tip_distances(tree)

    def brute(a: int, b: int) -> float:
        def up(x: int) -> dict[int, float]:
            out, d = {}, 0.0
            while x != -1:
                out[x] = d
                d += tree.length[x]
                x = tree.parent[x]
            return out

        ua, ub = up(a), up(b)
        return min(ua[n] + ub[n] for n in ua if n in ub)

    for node in range(tree.n_nodes):
        for k, tip in enumerate(tree.tips):
            assert dist[node, k] == pytest.approx(brute(node, tip))
        assert hi[node] - lo[node] >= 1


def test_decimal_year_and_augur_dates():
    assert decimal_year("2010-01-01") == pytest.approx(2010 + 0.5 / 365)
    assert decimal_year("2010-05") == pytest.approx(2010 + 4.5 / 12)
    assert decimal_year("2010") == 2010.5
    assert decimal_year("2010-06-XX") == decimal_year("2010-06")
    assert decimal_year("2010-XX-XX") == 2010.5
    assert decimal_year("XXXX-XX-XX") is None
    assert decimal_year("") is None
    assert augur_date("2010-05", "month") == "2010-05-XX"
    assert augur_date("2010", "year") == "2010-XX-XX"
    assert augur_date("2010-05-14", "day") == "2010-05-14"
    assert augur_date("", "unknown") == "?"


RATE = 0.002
YEARS = [2000.0, 2003.0, 2006.0, 2009.0, 2012.0, 2015.0]


def _clock_tree(years: list[float], *, odd_root: bool) -> str:
    lens = [RATE * (y - 1990.0) for y in years]
    names = [f"T{i}" for i in range(len(years))]
    if not odd_root:
        return "(" + ",".join(f"{n}:{x:.6f}" for n, x in zip(names, lens, strict=True)) + ");"
    rest = ",".join(f"{n}:{x:.6f}" for n, x in list(zip(names, lens, strict=True))[1:])
    return f"({names[0]}:{lens[0]:.6f},({rest}):0.0);"  # same star, root on a tip edge


@pytest.mark.parametrize("odd_root", [False, True])
def test_best_root_finds_the_clock_whatever_the_given_root(odd_root):
    tree = parse_newick(_clock_tree(YEARS, odd_root=odd_root))
    order = [int(tree.name[t][1:]) for t in tree.tips]
    dates = [YEARS[i] for i in order]
    used, d = best_root_to_tip(tree, dates)
    gate = evaluate_gate(used, d)
    assert gate.passed and gate.r_squared == pytest.approx(1.0, abs=1e-6)
    assert gate.slope == pytest.approx(RATE, rel=1e-3)


def test_shuffled_dates_fail_the_gate():
    tree = parse_newick(_clock_tree(YEARS, odd_root=False))
    shuffled = [2012.0, 2000.0, 2015.0, 2003.0, 2009.0, 2006.0]
    used, d = best_root_to_tip(tree, shuffled)
    assert not evaluate_gate(used, d, min_r2=0.9).passed


@pytest.fixture
def workspace(config_path, tmp_path):
    settings = load_settings(
        config_path,
        "dev",
        {
            "paths": {
                "raw": str(tmp_path / "raw"), "interim": str(tmp_path / "interim"),
                "processed": str(tmp_path / "processed"), "reference": str(tmp_path),
                "results": str(tmp_path / "results"), "auspice": str(tmp_path / "auspice"),
                "reports": str(tmp_path / "reports"),
            },
            "phylo": {"refine": {"enabled": True}},
        },
    )  # fmt: skip
    processed = tmp_path / "processed" / "kenya"
    write_fasta(processed / "aln.tree.fasta", [(f"T{i}", "ACGT") for i in range(6)], wrap=0)
    return settings, tmp_path


def _meta(tmp_path: Path, years: list[float]) -> None:
    write_tsv(
        tmp_path / "processed/kenya/metadata.final.tsv",
        [
            {"strain": f"T{i}", "date": f"{int(y)}-06", "date_precision": "month",
             "subtype_called": "A1"}
            for i, y in enumerate(years)
        ],
    )  # fmt: skip


def _fake_iqtree(tree_text: str):
    def runner(command):
        prefix = Path(command[command.index("--prefix") + 1])
        prefix.with_suffix(".treefile").write_text(tree_text + "\n")
        prefix.with_suffix(".log").write_text("log\n")

    return runner


def test_tree_command_and_outputs(workspace):
    settings, tmp = workspace
    cmd = iqtree_command(settings, tmp / "a.fasta", tmp / "p")
    assert cmd[0] == "iqtree2" and cmd[cmd.index("-m") + 1] == settings.phylo.tree.model
    assert ("-B" in cmd) == (settings.phylo.tree.bootstrap == "ufboot")
    assert cmd[cmd.index("-T") + 1] == str(settings.threads)
    out = run_tree(settings, repo_root=".", runner=_fake_iqtree(_clock_tree(YEARS, odd_root=False)))
    assert (tmp / "results/kenya/phylo/tree.nwk").exists()
    assert out["scope"] == "kenya"
    with pytest.raises(ContractError):
        run_tree(settings, repo_root=".", runner=lambda c: None)  # no treefile written


def _refine(workspace, years, tree_text):
    settings, tmp = workspace
    _meta(tmp, years)
    run_tree(settings, repo_root=".", runner=_fake_iqtree(tree_text))
    calls = []
    result = run_refine(settings, repo_root=".", runner=calls.append)
    gate = json.loads((tmp / "results/kenya/phylo/temporal_gate.json").read_text())
    return result, calls[0], gate, tmp


def test_refine_makes_a_time_tree_when_the_gate_passes(workspace):
    _, command, gate, tmp = _refine(workspace, YEARS, _clock_tree(YEARS, odd_root=False))
    assert gate["pooled"]["passed"] and gate["time_tree"] and "--timetree" in command
    assert "--keep-root" not in command
    assert (tmp / "results/kenya/phylo/metadata.augur.tsv").read_text().splitlines()[
        1
    ] == "T0\t2000-06-XX"


def test_refine_falls_back_to_a_divergence_tree_and_says_why(workspace):
    shuffled = [2012.0, 2000.0, 2015.0, 2003.0, 2009.0, 2006.0]
    tree_text = _clock_tree(YEARS, odd_root=False)  # clock follows YEARS, dates are shuffled
    settings, tmp = workspace
    strict = settings.model_copy(
        update={
            "phylo": settings.phylo.model_copy(
                update={
                    "temporal_gate": settings.phylo.temporal_gate.model_copy(
                        update={"min_temporal_r2": 0.9}
                    )
                }
            )
        }
    )
    _, command, gate, _ = _refine((strict, tmp), shuffled, tree_text)
    assert not gate["pooled"]["passed"] and not gate["time_tree"]
    assert "--keep-root" in command and "--timetree" not in command
    assert "R^2" in gate["note"]


def _pairwise(tree):
    d = tip_distances(tree)[0][tree.tips]
    names = [tree.name[t] for t in tree.tips]
    order = np.argsort(names)
    return d[np.ix_(order, order)]


def test_rerooting_keeps_distances_and_support_labels():
    tree = parse_newick("((A:1,B:2)95:0.5,(C:3,D:4)80:1.5,E:2);")
    for edge in [(2, 0.4), (1, 0.2), midpoint_edge(tree)]:
        text = to_newick(tree, edge)
        again = parse_newick(text)
        assert np.allclose(_pairwise(tree), _pairwise(again))
        assert "95" in text and "80" in text
    # midpoint of the longest path D-B (8) lies at the C/D clade node
    assert to_newick(tree, midpoint_edge(tree)).startswith("((C:3,D:4)80:0")


def test_permutation_p_separates_clock_from_noise():
    tree = parse_newick(_clock_tree(YEARS, odd_root=False))
    order = [int(tree.name[t][1:]) for t in tree.tips]
    dates = [YEARS[i] for i in order]
    used, d = best_root_to_tip(tree, dates)
    r2 = evaluate_gate(used, d).r_squared
    assert permutation_p(tree, dates, r2, n=99, seed=1) <= 0.05
    noisy = [2012.0, 2000.0, 2015.0, 2003.0, 2009.0, 2006.0]
    used, d = best_root_to_tip(tree, noisy)
    r2n = evaluate_gate(used, d).r_squared
    assert permutation_p(tree, noisy, r2n, n=99, seed=1) > 0.05


def test_failed_gate_roots_at_the_midpoint_and_reports_undated(workspace):
    shuffled = [2012.0, 2000.0, 2015.0, 2003.0, 2009.0, 2006.0]
    settings, tmp = workspace
    strict = settings.model_copy(
        update={
            "phylo": settings.phylo.model_copy(
                update={
                    "temporal_gate": settings.phylo.temporal_gate.model_copy(
                        update={"min_temporal_r2": 0.9}
                    )
                }
            )
        }
    )
    _, command, gate, _ = _refine((strict, tmp), shuffled, _clock_tree(YEARS, odd_root=False))
    assert command[command.index("--tree") + 1].endswith("tree.midpoint.nwk")
    assert (tmp / "results/kenya/phylo/tree.midpoint.nwk").exists()
    assert gate["n_tips"] == 6 and gate["n_undated"] == 0


def test_subtype_trees_get_their_own_inputs_gate_and_directory(workspace):
    settings, tmp = workspace
    settings = settings.model_copy(
        update={
            "phylo": settings.phylo.model_copy(
                update={
                    "refine": settings.phylo.refine.model_copy(update={"subtype_trees": ["A1"]})
                }
            )
        }
    )
    write_tsv(
        tmp / "processed/kenya/metadata.final.tsv",
        [
            {"strain": f"T{i}", "date": f"{int(y)}-06", "date_precision": "month",
             "subtype_called": "A1" if i < 5 else "D"}
            for i, y in enumerate(YEARS)
        ],
    )  # fmt: skip
    five = _clock_tree(YEARS[:5], odd_root=False)
    calls = []

    def runner(command):
        calls.append(command)
        if command[0] != "augur":
            _fake_iqtree(five)(command)

    from virotrace.phylo.stage import run_subtype_trees

    out = run_subtype_trees(settings, repo_root=".", runner=runner)
    assert out["subtype_trees"]["A1"]["time_tree"] is True
    assert (tmp / "processed/kenya/aln.tree.A1.fasta").read_text().count(">") == 5
    assert (tmp / "results/kenya/phylo/A1/temporal_gate.json").exists()
    assert "--timetree" in calls[-1]
    with pytest.raises(ContractError):  # D has one tip: no tree to build
        run_tree(settings, repo_root=".", runner=runner, subtype="D")
