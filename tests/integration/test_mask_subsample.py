"""Mask and subsample stages: gap-masked columns, and a tree set with a clean waterfall."""

from __future__ import annotations

import json

import pytest

from virotrace.curate.mask import masked_columns, run_mask
from virotrace.curate.subsample_stage import run_subsample
from virotrace.errors import ContractError
from virotrace.io.fasta import read_fasta, write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings


def test_masked_columns_uses_gap_and_n_only():
    seqs = ["AC-T", "AC-T", "ANGT", "AC-R"]
    # column 2 is gap in 3/4, column 1 is N in 1/4; R (ambiguity) is not missing
    assert masked_columns(seqs, 0.5).tolist() == [False, False, True, False]
    assert masked_columns(seqs, 0.9).tolist() == [False, False, False, False]


def test_masked_columns_rejects_ragged_input():
    with pytest.raises(ContractError):
        masked_columns(["AC", "ACG"], 0.5)


def _row(strain, isolate, year, precision, study="S1"):
    return {
        "strain": strain, "isolate": isolate, "year": year, "date": f"{year}-05" if year else "",
        "date_precision": precision, "country_code": "KEN", "ambiguity_frac": 0.0,
        "clone_derived": False, "study_key": study,
    }  # fmt: skip


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
            "subsample": {"max_per_stratum": 2, "year_bin_width": 5},
        },
    )  # fmt: skip
    rows = [
        _row("s1", "RV217_1_v01", 2010, "month"),
        _row("s2", "RV217_1_v02", 2011, "month"),  # same person as s1: dropped by policy
        _row("s3", "x3", 2010, "month"),
        _row("s4", "x4", 2011, "day"),
        _row("s5", "x5", 2012, "day"),  # third in the 2010-2014 A1 stratum: capped
        _row("s6", "x6", 2012, "year"),  # year-only: out under strict
        _row("s7", "x7", 2013, "month"),  # recombinant: out
        _row("s8", "x8", 2009, "month"),
    ]
    aln = [(r["strain"], "ACGT" + "A" * 4 + "-" * 2) for r in rows]
    interim = tmp_path / "interim" / "kenya"
    write_fasta(interim / "aln.pol.unmasked.fasta", aln, wrap=0)
    write_tsv(interim / "qc.passed.tsv", [{"strain": r["strain"]} for r in rows], ["strain"])
    write_tsv(interim / "metadata.valid.tsv", rows)
    write_tsv(
        interim / "subtype_calls.tsv",
        [{"strain": r["strain"], "subtype_called": "A1", "subtype_support": 1.0} for r in rows],
    )
    write_tsv(
        interim / "recombination.tsv",
        [
            {"strain": r["strain"], "recomb_status": "confirmed" if r["strain"] == "s7" else "none",
             "exclude_from_tree": r["strain"] == "s7"}
            for r in rows
        ],
    )  # fmt: skip
    return settings, tmp_path


def test_mask_then_subsample_gives_a_clean_waterfall(workspace):
    settings, tmp = workspace
    settings = settings.model_copy(
        update={"qc": settings.qc.model_copy(update={"individual_id_patterns": [r"^(RV217_\d+)_"]})}
    )
    masked = run_mask(settings, repo_root=".")
    assert masked["n_masked"] == 2 and masked["kept_ranges"] == [[1, 8]]
    assert {q[-2:] for _, q in read_fasta(tmp / "interim/kenya/aln.pol.masked.fasta")} == {"NN"}

    summary = run_subsample(settings, repo_root=".")
    wf = summary["waterfall"]
    assert wf["qc_passed"] == 8
    assert wf["after_individual_policy"] == 7  # s2 dropped
    assert wf["after_recombination_policy"] == 6  # s7 dropped
    assert wf["after_date_requirement"] == 5  # s6 (year-only) dropped
    assert wf["retained"] == 3  # 2010-2014 bin: 4 candidates capped at 2; 2005-2009 bin: s8
    final = read_tsv(tmp / "processed/kenya/metadata.final.tsv")
    assert len(final) == 3 and "s8" in {r["strain"] for r in final}
    assert "s2" not in {r["strain"] for r in final} and "s7" not in {r["strain"] for r in final}
    tree = {s for s, _ in read_fasta(tmp / "processed/kenya/aln.tree.fasta")}
    assert tree == {r["strain"] for r in final}
    on_disk = json.loads((tmp / "results/kenya/qc/subsample_summary.json").read_text())
    assert on_disk["year_span"][0] <= on_disk["year_span"][1]


def test_year_only_sensitivity_keeps_year_precision(workspace, config_path):
    settings, _ = workspace
    settings = settings.model_copy(
        update={"subsample": settings.subsample.model_copy(update={"date_requirement": "year"})}
    )
    run_mask(settings, repo_root=".")
    summary = run_subsample(settings, repo_root=".")
    assert (
        summary["waterfall"]["after_date_requirement"]
        == summary["waterfall"]["after_recombination_policy"]
    )
