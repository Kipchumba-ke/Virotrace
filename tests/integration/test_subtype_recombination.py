"""Subtype placement and the recombination screen on a synthetic panel with known truth."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pytest

from virotrace.curate.recombination import (
    classify,
    excluded_from_tree,
    run_recombination,
    segments_from_windows,
)
from virotrace.curate.subtype import (
    UNASSIGNED,
    agreement,
    build_panel,
    encode,
    panel_label,
    place,
    run_subtype,
)
from virotrace.errors import ContractError
from virotrace.io.fasta import write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings

POL_START, POL_END = 2085, 5096
LENGTH = POL_END - POL_START + 1
RNG = random.Random(11)
ROOT = "".join(RNG.choice("ACGT") for _ in range(LENGTH))


def mutate(seq: str, rate: float, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in seq)


# three "subtypes", each 12% from the root, two panel members each at 3% from their subtype
SUBTYPES = {name: mutate(ROOT, 0.12, i) for i, name in enumerate(["A1", "C", "D"])}
PANEL = [
    (f"{name}.KE.00.ref{k}.AB{i}{k}", mutate(seq, 0.03, 100 + 10 * i + k))
    for i, (name, seq) in enumerate(SUBTYPES.items())
    for k in range(2)
]


def cfg(config_path, **recomb):
    return load_settings(config_path, "dev", {"recombination": recomb} if recomb else None)


@pytest.fixture
def panel():
    return build_panel(PANEL)


def test_panel_labels_follow_lanl_headers():
    assert panel_label("A1.KE.00.name.AB123") == "A1"
    assert panel_label("01_AE.TH.90.CM240.U54771") == "CRF01_AE"
    with pytest.raises(ContractError):
        panel_label(".x")


def test_panel_needs_two_labels_and_equal_widths():
    with pytest.raises(ContractError):
        build_panel([("A1.x.1", "ACGT"), ("A1.x.2", "ACGT")])
    with pytest.raises(ContractError):
        build_panel([("A1.x.1", "ACGT"), ("C.x.1", "ACG")])


def test_a_query_is_placed_with_its_subtype_and_high_support(panel):
    rng = np.random.default_rng(1)
    for name, seq in SUBTYPES.items():
        placed = place(encode(mutate(seq, 0.04, 7)), panel, n_bootstrap=100, min_sites=300, rng=rng)
        assert placed.label == name and placed.support >= 0.95


def test_too_little_sequence_gives_no_placement(panel):
    short = "-" * 2500 + SUBTYPES["A1"][2500:2700] + "-" * (LENGTH - 2700)
    assert (
        place(encode(short), panel, n_bootstrap=10, min_sites=300, rng=np.random.default_rng(0))
        is None
    )


def test_ambiguous_codes_and_gaps_are_not_compared(panel):
    seq = list(SUBTYPES["C"])
    for i in range(0, 1500, 3):
        seq[i] = "N"
    placed = place(
        encode("".join(seq)), panel, n_bootstrap=20, min_sites=300, rng=np.random.default_rng(2)
    )
    assert placed.label == "C" and placed.n_sites < LENGTH - 400


def test_agreement_handles_group_labels_and_crfs():
    assert agreement("A1", "A") == "agree"
    assert agreement("A", "A1") == "agree"
    assert agreement("C", "A1") == "disagree"
    assert agreement("CRF01_AE", "CRF01_AE") == "agree"
    assert agreement("D", "") == "unreported"
    assert agreement(UNASSIGNED, "D") == "not_comparable"


def test_mosaic_is_suspected_and_confirmed_with_segments(panel, config_path):
    recomb = cfg(config_path).recombination
    mosaic = mutate(SUBTYPES["A1"][:1500] + SUBTYPES["D"][1500:], 0.03, 5)
    status, segments, _, _ = classify(encode(mosaic), panel, recomb, rng=np.random.default_rng(3))
    assert status == "confirmed"
    assert [s[0] for s in segments] == ["A1", "D"]
    assert 1300 <= segments[0][2] <= 1700  # the breakpoint sits near column 1500


def test_pure_sequences_are_not_flagged(panel, config_path):
    recomb = cfg(config_path).recombination
    for name, seq in SUBTYPES.items():
        status, segments, _, _ = classify(
            encode(mutate(seq, 0.04, 9)), panel, recomb, rng=np.random.default_rng(4)
        )
        assert status == "none" and [s[0] for s in segments] == [name]


def test_a_sequence_too_short_for_two_segments_is_not_testable(panel, config_path):
    recomb = cfg(config_path).recombination
    short = "-" * 1000 + SUBTYPES["C"][1000:1250] + "-" * (LENGTH - 1250)
    status, *_ = classify(encode(short), panel, recomb, rng=np.random.default_rng(5))
    assert status == "not_testable"


def test_segments_need_the_minimum_span_and_policies_are_independent():
    windows = [(0, 300, "A1"), (100, 400, "A1"), (200, 500, "D"), (300, 600, None)]
    assert segments_from_windows(windows, min_segment=400) == [("A1", 0, 400)]
    assert segments_from_windows(windows, min_segment=200) == [("A1", 0, 400), ("D", 200, 500)]
    assert excluded_from_tree("confirmed", "exclude_confirmed")
    assert not excluded_from_tree("suspected", "exclude_confirmed")
    assert excluded_from_tree("suspected", "exclude_suspected")
    assert not excluded_from_tree("confirmed", "keep_all")


# --- the stages ----------------------------------------------------------------------


def _identity_aligner():
    """Fake MAFFT for sequences already in reference coordinates (pad = same length)."""
    from virotrace.io.fasta import read_fasta

    def runner(ref_path: Path, frag_path: Path, threads: int) -> str:
        ((ref_name, ref_seq),) = read_fasta(ref_path)
        out = [f">{ref_name}\n{ref_seq}\n"]
        out += [f">{name}\n{seq}\n" for name, seq in read_fasta(frag_path)]
        return "".join(out)

    return runner


def _hxb2_genbank(path: Path) -> None:
    seq = "A" * (POL_START - 1) + ROOT + "C" * (9719 - POL_END)
    lines = ["LOCUS       K03455                9719 bp    DNA     linear   VRL 01-JAN-2000"]
    lines += ["VERSION     K03455.1", "FEATURES             Location/Qualifiers"]
    lines += ["     source          1..9719", "ORIGIN"]
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60].lower()
        lines.append(f"{i + 1:>9} " + " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10)))
    lines.append("//")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def workspace(config_path, tmp_path):
    _hxb2_genbank(tmp_path / "hxb2.gb")
    write_fasta(tmp_path / "panel.fasta", PANEL, wrap=0)
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
            "subtype": {"reference_panel": "panel.fasta"},
        },
    )  # fmt: skip
    queries = [
        ("PURE_A1", mutate(SUBTYPES["A1"], 0.04, 21)),
        ("PURE_C", mutate(SUBTYPES["C"], 0.04, 22)),
        ("MOSAIC", mutate(SUBTYPES["A1"][:1500] + SUBTYPES["D"][1500:], 0.03, 23)),
        ("SHORT", "-" * 1000 + SUBTYPES["C"][1000:1150] + "-" * (LENGTH - 1150)),
        ("FILTERED", mutate(SUBTYPES["D"], 0.04, 24)),
    ]
    interim = tmp_path / "interim" / "kenya"
    write_fasta(interim / "aln.pol.unmasked.fasta", queries, wrap=0)
    write_tsv(interim / "qc.passed.tsv", [{"strain": s} for s, _ in queries[:4]], ["strain"])
    write_tsv(
        interim / "metadata.valid.tsv",
        [{"strain": s, "subtype_reported": r} for s, r in
         [("PURE_A1", "A"), ("PURE_C", "D"), ("MOSAIC", ""), ("SHORT", "C"), ("FILTERED", "D")]],
        ["strain", "subtype_reported"],
    )  # fmt: skip
    return settings, tmp_path


def test_subtype_stage_calls_reports_agreement_and_skips_unpassed(workspace):
    settings, tmp = workspace
    summary = run_subtype(settings, runner=_identity_aligner())
    calls = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "subtype_calls.tsv")}
    assert set(calls) == {"PURE_A1", "PURE_C", "MOSAIC", "SHORT"}  # FILTERED did not pass QC
    assert (
        calls["PURE_A1"]["subtype_called"] == "A1"
        and calls["PURE_A1"]["subtype_agreement"] == "agree"
    )
    assert (
        calls["PURE_C"]["subtype_called"] == "C"
        and calls["PURE_C"]["subtype_agreement"] == "disagree"
    )
    assert calls["SHORT"]["subtype_called"] == UNASSIGNED
    assert summary["n_sequences"] == 4 and summary["panel"]["n_sequences"] == 6
    assert summary["disagreement_pairs"] == {"C|D": 1}
    assert summary["unassigned_by_n_sites"] == {"under_600": 1, "600_or_more": 0}


def test_recombination_stage_flags_the_mosaic_and_applies_the_tree_policy(workspace):
    settings, tmp = workspace
    run_subtype(settings, runner=_identity_aligner())
    summary = run_recombination(settings)
    flags = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "recombination.tsv")}
    assert (
        flags["MOSAIC"]["recomb_status"] == "confirmed"
        and flags["MOSAIC"]["exclude_from_tree"] == "true"
    )
    assert (
        flags["PURE_A1"]["recomb_status"] == "none"
        and flags["PURE_A1"]["exclude_from_tree"] == "false"
    )
    assert flags["SHORT"]["recomb_status"] == "not_testable"
    assert summary["segment_patterns"] == {"A>D": 1}
    assert summary["n_excluded_from_tree"] == 1 and summary["drm_policy"] == "keep_all"


def test_missing_panel_says_where_to_get_it(config_path, workspace):
    settings, tmp = workspace
    (tmp / "panel.fasta").unlink()
    with pytest.raises(ContractError, match="LANL"):
        run_subtype(settings, runner=_identity_aligner())


def test_recombination_requires_the_subtype_stage(workspace):
    settings, _ = workspace
    with pytest.raises(ContractError, match="subtype"):
        run_recombination(settings)


def test_ref_prefixed_lanl_headers_and_panel_filters(config_path):
    from virotrace.curate.subtype import select_panel

    assert panel_label("Ref.A1.AU.03.PS1044_Day0.DQ676872") == "A1"
    assert panel_label("Ref.01_AE.TH.90.CM240.U54771") == "CRF01_AE"
    rows = [(h, "ACGT") for h in (
        "Ref.A1.AU.03.a.X1", "Ref.CPZ.CM.01.b.X2", "Ref.O.CM.01.c.X3", "Ref.GOR.CM.01.d.X4",
        "Ref.01_AE.TH.90.e.X5", "Ref.02_AG.NG.x.f.X6", "Ref.94_cpx.CD.1.g.X7",
    )]  # fmt: skip
    cfg = load_settings(config_path, "dev").subtype
    assert [n.split(".")[1] for n, _ in select_panel(rows, cfg)] == [
        "A1",
        "01_AE",
        "02_AG",
        "94_cpx",
    ]
    narrowed = load_settings(config_path, "dev", {"subtype": {"include_crfs": ["01_AE"]}}).subtype
    assert [n.split(".")[1] for n, _ in select_panel(rows, narrowed)] == ["A1", "01_AE"]


def test_panel_self_check_recovers_labels_and_reports_failures(panel):
    from virotrace.curate.subtype import panel_self_check

    result = panel_self_check(
        panel, n_bootstrap=50, min_sites=300, min_support=0.7, rng=np.random.default_rng(0)
    )
    assert result["n_tested"] == 6 and result["n_recovered"] == 6 and result["failures"] == []
    # a mislabelled entry is a finding, not a crash
    rows = [(n.replace("D.KE", "C.KE") if n.startswith("D.KE.00.ref0") else n, s) for n, s in PANEL]
    bad = panel_self_check(
        build_panel(rows),
        n_bootstrap=50,
        min_sites=300,
        min_support=0.7,
        rng=np.random.default_rng(0),
    )
    assert bad["n_failures"] >= 1 and bad["failures"][0]["expected"] == "C"


def test_stage_summary_carries_the_self_check(workspace):
    settings, _ = workspace
    summary = run_subtype(settings, runner=_identity_aligner())
    assert summary["panel_self_check"]["recovery_rate"] == 1.0


# --- CRFs near a pure subtype must not split the call ---------------------------------

CRF_NEAR_A1 = [
    (f"50_A1D.KE.10.crf{k}.CC{k}", mutate(SUBTYPES["A1"], 0.015, 300 + k)) for k in range(3)
]


@pytest.fixture
def crf_workspace(workspace):
    settings, tmp = workspace
    write_fasta(tmp / "panel.fasta", PANEL + CRF_NEAR_A1, wrap=0)
    return settings, tmp


def test_crfs_near_a_pure_subtype_do_not_leave_it_unassigned(crf_workspace):
    settings, tmp = crf_workspace
    summary = run_subtype(settings, runner=_identity_aligner())
    calls = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "subtype_calls.tsv")}
    assert calls["PURE_A1"]["subtype_called"] == "A1"  # not split with CRF50_A1D
    assert float(calls["PURE_A1"]["subtype_support"]) >= 0.9
    assert summary["panel"]["n_pure"] == 6 and summary["panel"]["n_crf_labels"] == 1
    assert "panel_self_check_with_crfs" in summary


def test_crf_label_is_annotation_only_and_can_be_switched_off(crf_workspace, config_path):
    settings, tmp = crf_workspace
    run_subtype(settings, runner=_identity_aligner())
    calls = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "subtype_calls.tsv")}
    assert calls["PURE_C"]["crf_label"] == ""
    off = load_settings(
        config_path,
        "dev",
        {
            "paths": settings.paths,
            "subtype": {"reference_panel": "panel.fasta", "annotate_crfs": False},
        },
    )
    run_subtype(off, runner=_identity_aligner())
    again = read_tsv(tmp / "interim" / "kenya" / "subtype_calls.tsv")
    assert all(r["crf_label"] == "" for r in again)


def test_recombination_parents_are_pure_groups_even_with_crfs_in_the_panel(crf_workspace):
    settings, tmp = crf_workspace
    run_subtype(settings, runner=_identity_aligner())
    run_recombination(settings)
    flags = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "recombination.tsv")}
    assert flags["MOSAIC"]["recomb_status"] == "confirmed"
    assert flags["PURE_A1"]["recomb_status"] == "none"  # CRF50_A1D is not a competing parent
