"""TN93 distances (against a scalar reference) and the cluster stage."""

from __future__ import annotations

import json
import math
import random
import shutil

import pytest

from virotrace.io.fasta import write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.phylo.clusters import mantel_haenszel, run_clusters
from virotrace.phylo.tn93 import components, tn93_edges
from virotrace.settings import load_settings


def tn93_reference(a: str, b: str) -> float:
    """Straight from the formula, one pair at a time."""
    pairs = [(x, y) for x, y in zip(a, b, strict=True) if x in "ACGT" and y in "ACGT"]
    n = len(pairs)
    freq = {k: 0.0 for k in "ACGT"}
    for x, y in pairs:
        freq[x] += 1
        freq[y] += 1
    pa, pc, pg, pt = (freq[k] / (2 * n) for k in "ACGT")
    pr, py = pa + pg, pc + pt
    purines, pyr = set("AG"), set("CT")
    p1 = sum(1 for x, y in pairs if x != y and {x, y} <= purines) / n
    p2 = sum(1 for x, y in pairs if x != y and {x, y} <= pyr) / n
    q = sum(1 for x, y in pairs if x != y and ({x, y} - purines and {x, y} - pyr)) / n
    t1 = 1 - pr * p1 / (2 * pa * pg) - q / (2 * pr)
    t2 = 1 - py * p2 / (2 * pc * pt) - q / (2 * py)
    t3 = 1 - q / (2 * pr * py)
    return (
        -2 * pa * pg / pr * math.log(t1)
        - 2 * pc * pt / py * math.log(t2)
        - 2 * (pr * py - pa * pg * py / pr - pc * pt * pr / py) * math.log(t3)
    )


def _mutate(seq: str, rng: random.Random, n: int) -> str:
    out = list(seq)
    for i in rng.sample(range(len(out)), n):
        out[i] = rng.choice([c for c in "ACGT" if c != out[i]])
    return "".join(out)


def test_tn93_matches_the_scalar_formula_including_ambiguity_and_gaps():
    rng = random.Random(5)
    base = "".join(rng.choice("ACGT") for _ in range(600))
    seqs = [base, _mutate(base, rng, 12), _mutate(base, rng, 40), _mutate(base, rng, 90)]
    seqs[1] = "N" * 20 + seqs[1][20:60] + "-" * 15 + seqs[1][75:]
    edges, compared = tn93_edges(seqs, threshold=1.0, min_overlap=300)
    assert compared == 6
    got = {(i, j): d for i, j, d in edges}
    for (i, j), d in got.items():
        assert d == pytest.approx(tn93_reference(seqs[i], seqs[j]), rel=1e-6)
    assert len(got) == 6


def test_identical_sequences_have_zero_distance_and_short_overlap_is_unknown():
    seq = "ACGT" * 100
    edges, compared = tn93_edges([seq, seq, "N" * 380 + "ACGT" * 5], 0.015, min_overlap=300)
    assert edges == [(0, 1, 0.0)] and compared == 1  # the short-overlap pairs are not compared


def test_components_are_single_linkage_and_sorted_largest_first():
    groups = components(6, [(0, 1, 0.0), (1, 2, 0.0), (4, 5, 0.0)])
    assert groups == [[0, 1, 2], [4, 5], [3]]


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
            "phylo": {"clusters": {"enabled": True, "threshold": 0.015, "min_overlap": 300}},
        },
    )  # fmt: skip
    shutil.copy("config/reference/sdrm_who.tsv", tmp_path / "sdrm_who.tsv")
    rng = random.Random(9)
    base = "".join(rng.choice("ACGT") for _ in range(600))
    other = "".join(rng.choice("ACGT") for _ in range(600))
    seqs = {
        "a1": base, "a2": _mutate(base, rng, 3), "a3": _mutate(base, rng, 4),
        "b1": other, "b2": _mutate(other, rng, 2),
        "c1": "".join(rng.choice("ACGT") for _ in range(600)),
    }  # fmt: skip
    processed = tmp_path / "processed" / "kenya"
    write_fasta(processed / "aln.tree.fasta", list(seqs.items()), wrap=0)
    write_tsv(
        processed / "metadata.final.tsv",
        [
            {"strain": s, "subtype_called": "A1", "study_key": "S1" if s[0] != "b" else "S2"}
            for s in seqs
        ],
    )
    drm = [
        {"strain": s, "gene": "RT", "position": 103, "ref_aa": "K", "residues": "N",
         "sdrm_residues": "N"}
        for s in ("a1", "a2", "c1")
    ]  # fmt: skip
    write_tsv(
        tmp_path / "results/kenya/drm/drm_table.tsv",
        drm,
        ["strain", "gene", "position", "ref_aa", "residues", "sdrm_residues"],
    )
    return settings, tmp_path


def test_cluster_stage_reports_clusters_and_sdrm_by_membership(workspace):
    settings, tmp = workspace
    summary = run_clusters(settings, repo_root=".")
    assert summary["n_clusters"] == 2 and summary["size_distribution"] == {2: 1, 3: 1}
    assert summary["n_clustered"] == 5 and summary["largest_cluster"] == 3
    assert summary["clusters_sharing_an_sdrm"] == 1  # a1 and a2 share RT:K103N
    assert summary["sdrm"]["clustered"]["with_sdrm"] == 2
    assert summary["sdrm"]["unclustered"] == {
        "n": 1, "with_sdrm": 1, "proportion": 1.0, "ci_low": 0.2065, "ci_high": 1.0,
    }  # fmt: skip
    rows = {r["strain"]: r for r in read_tsv(tmp / "results/kenya/phylo/clusters.tsv")}
    assert rows["a1"]["cluster_id"] == rows["a3"]["cluster_id"] == "CL001"
    assert rows["c1"]["cluster_id"] == "" and rows["c1"]["cluster_size"] == "1"
    bins = summary["edges_by_distance"]
    assert sum(b["n_pairs"] for b in bins.values()) == 4  # a1-a2, a1-a3, a2-a3, b1-b2
    assert all(b["cross_study"] == 0 for b in bins.values())
    assert summary["sdrm_by_treatment_context"]["unannotated"]["clustered"]["n"] == 5
    assert summary["near_identical_pairs"]["n"] >= 0
    assert (tmp / "results/kenya/qc/near_identical_pairs.tsv").exists()
    masked = summary["sensitivity_sdrm_masked"]
    assert masked["n_masked_columns"] > 0 and masked["n_clusters"] == 2
    on_disk = json.loads((tmp / "results/kenya/phylo/clusters_summary.json").read_text())
    assert on_disk["multi_study_clusters"] == 0


def test_mantel_haenszel_single_stratum_matches_woolf():
    out = mantel_haenszel([(10, 20, 5, 40)])
    se = math.sqrt(1 / 10 + 1 / 20 + 1 / 5 + 1 / 40)
    assert out["odds_ratio"] == 4.0
    assert out["ci_low"] == pytest.approx(math.exp(math.log(4) - 1.96 * se), abs=1e-3)
    assert out["ci_high"] == pytest.approx(math.exp(math.log(4) + 1.96 * se), abs=1e-3)
    assert mantel_haenszel([(0, 0, 0, 0)])["odds_ratio"] is None


def test_stratification_removes_a_confounded_association():
    # within each cohort clustered and unclustered have the same SDRM rate: OR must be 1
    assert mantel_haenszel([(9, 1, 18, 2), (1, 9, 2, 18)])["odds_ratio"] == 1.0
