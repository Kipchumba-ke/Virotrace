"""The drm-prevalence stage on a synthetic population with known answers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from virotrace.drm.classes import NNRTI_POSITIONS, NRTI_POSITIONS, drug_class_of
from virotrace.drm.prevalence_stage import run_prevalence, year_bin
from virotrace.drm.sdrm import load_sdrm_list
from virotrace.errors import ContractError
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings

REPO = Path(__file__).resolve().parents[2]
COV = [
    "strain", "gene", "coverage", "usable", "sdrm_position_coverage",
    "sdrm_position_coverage_sensitivity",
]  # fmt: skip
TABLE = [
    "strain", "gene", "position", "is_sdrm", "is_sdrm_sensitivity",
]  # fmt: skip


@pytest.fixture
def workspace(config_path, tmp_path):
    (tmp_path / "head.tsv").write_text(
        "gene\tposition\tresidues\nRT\t184\tV\nRT\t103\tN\nPR\t46\tI\n"
    )
    (tmp_path / "sens.tsv").write_text(
        "gene\tposition\tresidues\nRT\t184\tV\nRT\t103\tN\nPR\t46\tI\nRT\t318\tF\nIN\t148\tH\n"
    )
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
            "drm": {"sdrm": {
                "list": "head.tsv", "list_version": "head",
                "sensitivity_list": "sens.tsv", "sensitivity_list_version": "sens",
            }},
        },
    )  # fmt: skip
    strains = [f"S{i}" for i in range(48)]
    years = {s: 2001 + 5 * (i // 12) for i, s in enumerate(strains)}
    cov, table, meta = [], [], []
    for s in strains:
        for gene in ("PR", "RT", "IN"):
            usable = "false" if (s, gene) == ("S6", "RT") else "true"
            head = "" if gene == "IN" else "1.0"
            sens = "0.5" if (s, gene) == ("S5", "RT") else "1.0"
            cov.append(
                {"strain": s, "gene": gene, "coverage": "1.0", "usable": usable,
                 "sdrm_position_coverage": head, "sdrm_position_coverage_sensitivity": sens}
            )  # fmt: skip
        meta.append(
            {"strain": s, "year": years[s], "subtype_reported": "A1" if int(s[1:]) % 2 else "",
             "clone_derived": "true" if s == "S40" else "false"}
        )  # fmt: skip

    def hit(s, pos, head=True):
        table.append(
            {"strain": s, "gene": "RT", "position": pos, "is_sdrm": str(head).lower(),
             "is_sdrm_sensitivity": "true"}
        )  # fmt: skip

    for i in range(36, 48):
        hit(f"S{i}", 184)  # NRTI in the last bin (S40 is a clone, S41 fails QC)
    for i in range(24, 28):
        hit(f"S{i}", 103)  # NNRTI
    hit("S6", 184)  # RT unusable: must not count
    hit("S1", 318, head=False)  # only on the sensitivity list

    interim = tmp_path / "interim" / "kenya"
    write_tsv(
        interim / "metadata.valid.tsv",
        meta,
        ["strain", "year", "subtype_reported", "clone_derived"],
    )
    write_tsv(interim / "qc.passed.tsv", [{"strain": s} for s in strains if s != "S41"], ["strain"])
    drm = tmp_path / "results" / "kenya" / "drm"
    write_tsv(drm / "gene_coverage.tsv", cov, COV)
    write_tsv(drm / "drm_table.tsv", table, TABLE)
    return settings, tmp_path


def _rows(tmp, **where):
    rows = read_tsv(tmp / "results" / "kenya" / "drm" / "prevalence.tsv")
    return [r for r in rows if all(r[k] == v for k, v in where.items())]


def test_denominators_are_per_class_and_exclude_unusable_clone_and_failed_qc(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    nrti = _rows(tmp, list_role="headline", drug_class="NRTI", stratum_kind="overall")[0]
    # 48 - S6 (RT unusable) - S40 (clone) - S41 (failed QC)
    assert (nrti["numerator"], nrti["denominator"]) == ("10", "45")
    nnrti = _rows(tmp, list_role="headline", drug_class="NNRTI", stratum_kind="overall")[0]
    assert (nnrti["numerator"], nnrti["denominator"]) == ("4", "45")
    pi = _rows(tmp, list_role="headline", drug_class="PI", stratum_kind="overall")[0]
    assert (pi["numerator"], pi["denominator"]) == ("0", "46")  # S6's protease is fine


def test_a_class_without_list_positions_is_not_estimated(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    (insti,) = _rows(tmp, list_role="headline", drug_class="INSTI")
    assert insti["status"] == "no_list_positions" and insti["proportion"] == ""
    sens = _rows(tmp, list_role="sensitivity", drug_class="INSTI", stratum_kind="overall")[0]
    assert sens["status"] == "estimated" and sens["denominator"] == "46"


def test_sensitivity_list_has_its_own_denominators_and_hits(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    sens = _rows(tmp, list_role="sensitivity", drug_class="NNRTI", stratum_kind="overall")[0]
    # S5 does not cover the extra positions; S1's 318 counts only here
    assert (sens["numerator"], sens["denominator"]) == ("5", "44")


def test_year_bins_wilson_and_suppression(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    last = _rows(
        tmp, list_role="headline", drug_class="NRTI", stratum_kind="year_bin", stratum="2015-2019"
    )[0]
    assert (last["numerator"], last["denominator"]) == ("10", "10")
    assert float(last["ci_low"]) < 1.0 == float(last["ci_high"])
    subtype = _rows(tmp, list_role="headline", drug_class="PI", stratum_kind="subtype_reported")
    assert {r["stratum"] for r in subtype} == {"A1", "unassigned"}


def test_small_strata_are_suppressed_not_dropped(config_path, workspace):
    settings, tmp = workspace
    raised = load_settings(config_path, "dev", {"paths": settings.paths, "drm": settings.drm.model_dump(),
                                                 "prevalence": {"min_stratum_n": 30}})  # fmt: skip
    run_prevalence(raised)
    row = _rows(
        tmp, list_role="headline", drug_class="NRTI", stratum_kind="year_bin", stratum="2015-2019"
    )[0]
    assert row["status"] == "n_too_small" and row["proportion"] == ""


def test_trend_is_tested_and_labelled_with_its_caveat(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    trends = {
        (r["list_role"], r["drug_class"]): r
        for r in read_tsv(tmp / "results" / "kenya" / "drm" / "trend_tests.tsv")
    }
    nrti = trends[("headline", "NRTI")]
    assert nrti["n_bins"] == "4" and float(nrti["p_two_sided"]) < 0.001
    assert "sampling" in nrti["note"]
    assert ("headline", "INSTI") not in trends  # nothing to test


def test_summary_is_written_and_missing_inputs_name_the_stage(workspace, config_path, tmp_path):
    settings, tmp = workspace
    summary = run_prevalence(settings)
    assert summary["n_population"] == 46
    assert (
        json.loads((tmp / "results" / "kenya" / "drm" / "prevalence_summary.json").read_text())[
            "scope"
        ]
        == "kenya"
    )
    empty = load_settings(
        config_path, "dev", {"paths": {**settings.paths, "results": str(tmp_path / "nowhere")}}
    )
    with pytest.raises(ContractError, match="drm call"):
        run_prevalence(empty)


def test_drug_classes_cover_both_real_lists_and_year_bins():
    for name in ("sdrm_who.tsv", "sdrm_stanford_2026-09-25.tsv"):
        for gene, pos in load_sdrm_list(REPO / "config" / "reference" / name, "x").entries:
            assert drug_class_of(gene, pos)
    assert not NRTI_POSITIONS & NNRTI_POSITIONS
    assert year_bin(2012, 5) == "2010-2014"
    with pytest.raises(ContractError):
        drug_class_of("RT", 999)


def test_same_denominator_sensitivity_matches_headline_people(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    head = _rows(tmp, list_role="headline", drug_class="NNRTI", stratum_kind="overall")[0]
    same = _rows(
        tmp,
        list_role="sensitivity_headline_denominator",
        drug_class="NNRTI",
        stratum_kind="overall",
    )[0]
    assert same["denominator"] == head["denominator"] == "45"
    assert same["numerator"] == "5"  # the 4 headline hits plus S1's R318F-only hit
    (insti,) = _rows(tmp, list_role="sensitivity_headline_denominator", drug_class="INSTI")
    assert insti["status"] == "no_list_positions"


def test_trend_p_value_does_not_underflow_to_zero():
    from virotrace.drm.prevalence import cochran_armitage_trend

    _, p = cochran_armitage_trend([0, 50, 150, 300], [1000, 1000, 1000, 1000])
    assert 0 < p < 1e-30


def test_study_summary_names_the_dominant_studies(workspace):
    settings, tmp = workspace
    run_prevalence(settings)
    rows = read_tsv(tmp / "results" / "kenya" / "drm" / "study_summary.tsv")
    assert rows[0]["study_key"] == "unknown_study" and rows[0]["n_sequences"] == "46"
    assert rows[0]["NRTI"] == "10/45"


def test_one_per_individual_shrinks_the_population_and_annotations_override_the_hint(
    config_path, workspace
):

    settings, tmp = workspace
    meta_path = tmp / "interim" / "kenya" / "metadata.valid.tsv"
    rows = read_tsv(meta_path)
    for i, row in enumerate(rows):
        row["isolate"] = f"RV217_{i // 6}_v01_GP_{i}"  # eight "participants" of six sequences
        row["study_key"] = "REF:Test|study"
        row["date"] = f"2010-{(i % 6) + 1:02d}"
        row["ambiguity_frac"] = "0.0"
    write_tsv(meta_path, rows, list(rows[0]))
    (tmp / "annotations.tsv").write_text(
        "study_key\ttreatment_context\tnote\nREF:Test|study\texperienced\ttest\n"
    )
    custom = load_settings(
        config_path,
        "dev",
        {
            "paths": settings.paths,
            "drm": settings.drm.model_dump(),
            "prevalence": {"study_annotations": "annotations.tsv"},
        },
    )
    summary = run_prevalence(custom)
    assert summary["n_population"] == 8  # 48 - S41 (QC) leaves 47; one per participant = 8
    assert summary["n_removed_by_population_policy"]["individual_policy"] > 0
    statuses = {
        r["stratum"]
        for r in read_tsv(tmp / "results" / "kenya" / "drm" / "prevalence.tsv")
        if r["stratum_kind"] == "treatment_status"
    }
    assert statuses == {"experienced"}


def test_annotation_file_is_validated(tmp_path):
    from virotrace.curate.study import load_study_annotations

    assert load_study_annotations(tmp_path / "none.tsv") == {}
    bad = tmp_path / "bad.tsv"
    bad.write_text("study_key\ttreatment_context\nX\tsometimes\n")
    with pytest.raises(ContractError):
        load_study_annotations(bad)
    real = load_study_annotations(REPO / "config" / "reference" / "study_annotations.tsv")
    assert (
        real["REF:Rossenkhan|bayesian-estimation-of-hiv-acquisition-d"]["treatment_context"]
        == "naive"
    )


def test_bad_individual_pattern_is_rejected(config_path):
    from virotrace.errors import ConfigError

    with pytest.raises((ConfigError, ValueError)):
        load_settings(config_path, "dev", {"qc": {"individual_id_patterns": ["no group"]}})


def test_keyword_hint_never_stratifies_only_annotations_do(config_path, workspace):
    settings, tmp = workspace
    meta_path = tmp / "interim" / "kenya" / "metadata.valid.tsv"
    rows = read_tsv(meta_path)
    for row in rows:
        row["treatment_status"] = "naive"  # what the keyword hint would say
        row["study_key"] = "REF:Unreviewed|study"
    write_tsv(meta_path, rows, list(rows[0]))
    run_prevalence(settings)
    statuses = {
        r["stratum"]
        for r in read_tsv(tmp / "results" / "kenya" / "drm" / "prevalence.tsv")
        if r["stratum_kind"] == "treatment_status"
    }
    assert statuses == {"unannotated"}
    studies = read_tsv(tmp / "results" / "kenya" / "drm" / "study_summary.tsv")
    assert studies[0]["treatment_hint"] == "naive" and studies[0]["annotation"] == "none"


def test_prevalence_stratifies_by_the_called_subtype_when_calls_exist(workspace):
    settings, tmp = workspace
    strains = [r["strain"] for r in read_tsv(tmp / "interim" / "kenya" / "metadata.valid.tsv")]
    write_tsv(
        tmp / "interim" / "kenya" / "subtype_calls.tsv",
        [{"strain": s, "subtype_called": "D" if int(s[1:]) % 3 else "A1"} for s in strains],
        ["strain", "subtype_called"],
    )
    run_prevalence(settings)
    kinds = {
        r["stratum_kind"] for r in read_tsv(tmp / "results" / "kenya" / "drm" / "prevalence.tsv")
    }
    assert "subtype_called" in kinds and "subtype_reported" not in kinds
