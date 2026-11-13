"""The report: laid out from stage outputs, limitations always present."""

from __future__ import annotations

import json
import shutil

import pytest

from virotrace.errors import ContractError
from virotrace.export.report import limitations, run_report
from virotrace.io.tables import write_tsv
from virotrace.settings import load_settings

PREV = [
    "list_role", "list_version", "stratum_kind", "stratum", "drug_class", "numerator",
    "denominator", "proportion", "ci_low", "ci_high", "status",
]  # fmt: skip


def _row(kind, stratum, klass, num, den, status="estimated", role="headline"):
    ok = status == "estimated"
    return {
        "list_role": role, "list_version": "2009", "stratum_kind": kind, "stratum": stratum,
        "drug_class": klass, "numerator": num, "denominator": den,
        "proportion": 0.25 if ok else "", "ci_low": 0.1 if ok else "",
        "ci_high": 0.4 if ok else "", "status": status,
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
            "phylo": {"refine": {"subtype_trees": ["A1"]}},
        },
    )  # fmt: skip
    results = tmp_path / "results/kenya"
    write_tsv(
        results / "drm/prevalence.tsv",
        [
            _row("overall", "all", "NRTI", 25, 100),
            _row("overall", "all", "NNRTI", 3, 5, "n_too_small"),
            _row("subtype_called", "A1", "NRTI", 20, 80),
            _row("overall", "all", "NRTI", 30, 100, role="sensitivity"),
        ],
        PREV,
    )
    write_tsv(
        results / "drm/trend_tests.tsv",
        [{"list_role": "headline", "list_version": "2009", "drug_class": "NRTI", "n_bins": 3,
          "first_bin": "2008-2012", "last_bin": "2018-2022", "z": 1.2, "p_two_sided": 0.2,
          "note": ""}],
        ["list_role", "list_version", "drug_class", "n_bins", "first_bin", "last_bin", "z",
         "p_two_sided", "note"],
    )  # fmt: skip
    (results / "qc").mkdir(parents=True)
    (results / "qc/qc_summary.json").write_text(
        json.dumps({"n_passed": 90, "waterfall": {"QC-0": 10, "QC-2": 3},
                    "clone_policy": "exclude", "individual_policy": "one_per_individual",
                    "n_removed_by_population_policy": {"clone_policy": 5},
                    "corpus_after_qc_and_clone_policy": {
                        "n_total": 80, "n_time_tree_eligible_strict": 60,
                        "n_time_tree_eligible_with_year_only": 70}})
    )  # fmt: skip
    (results / "qc/subsample_summary.json").write_text(
        json.dumps({"waterfall": {"retained": 60}, "year_span": [2008, 2020],
                    "by_subtype": {"A1": 40, "D": 20}, "by_year": {"2010": 30, "2015": 30}})
    )  # fmt: skip
    gate = {"n": 60, "slope": 0.001, "r_squared": 0.05, "permutation_p": 0.3, "passed": False}
    phylo = results / "phylo"
    (phylo / "A1").mkdir(parents=True)
    (phylo / "temporal_gate.json").write_text(
        json.dumps({"pooled": gate, "time_tree": False, "by_subtype": {"A1": gate}})
    )
    (phylo / "A1/temporal_gate.json").write_text(
        json.dumps({"pooled": {**gate, "passed": True, "r_squared": 0.12}, "time_tree": True})
    )
    (phylo / "clusters_summary.json").write_text(
        json.dumps({
            "threshold": 0.015, "n_clustered": 10, "n_sequences": 60,
            "proportion_clustered": 0.1667, "n_clusters": 4, "largest_cluster": 3,
            "sdrm": {"clustered": {"n": 10, "with_sdrm": 6, "proportion": 0.6},
                     "unclustered": {"n": 50, "with_sdrm": 20, "proportion": 0.4}},
            "sdrm_cluster_odds_ratio_adjusted_for_treatment_context":
                {"odds_ratio": 1.3, "ci_low": 0.8, "ci_high": 2.2},
        })
    )  # fmt: skip
    (tmp_path / "docs").mkdir()
    shutil.copy("docs/04-bioinformatics-methodology.md", tmp_path / "docs/")
    return settings, tmp_path


def test_report_has_headline_estimates_suppression_gates_and_clusters(workspace):
    settings, tmp = workspace
    summary = run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert summary["report"].endswith("reports/kenya/index.html")
    assert "25.0% (10.0% to 40.0%), 25/100" in html
    assert "suppressed (n=5)" in html
    assert "Stanford" in html and "By subtype called" in html
    assert "divergence tree" in html and "A1 tree" in html and "time tree" in html
    assert "1.3 (95% interval 0.8 to 2.2)" in html
    assert "Not a clinical tool" in html


def test_limitations_come_from_the_methods_doc_and_are_all_in_the_report(workspace):
    settings, tmp = workspace
    items = limitations(tmp)
    assert len(items) >= 14 and items[0].startswith("Convenience sampling")
    run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert all(item.split(" ")[0] in html for item in items)
    limits = html.split("<h2>Limitations</h2>")[1].split("<h2>References</h2>")[0]
    assert limits.count("<li>") == len(items)


def test_report_without_phylo_outputs_says_so_and_missing_prevalence_names_the_stage(workspace):
    settings, tmp = workspace
    (tmp / "results/kenya/phylo/temporal_gate.json").unlink()
    run_report(settings, repo_root=tmp)
    assert "Not generated for this build." in (tmp / "reports/kenya/index.html").read_text()
    (tmp / "results/kenya/drm/prevalence.tsv").unlink()
    with pytest.raises(ContractError, match="drm prevalence"):
        run_report(settings, repo_root=tmp)


def test_report_html_escapes_study_text(workspace):
    settings, tmp = workspace
    write_tsv(
        tmp / "results/kenya/drm/prevalence.tsv",
        [_row("study", "<script>x</script>", "NRTI", 25, 100)],
        PREV,
    )
    run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert "<script>x" not in html and "&lt;script&gt;" in html


def test_report_shows_clock_sensitivity_and_the_tree_based_subtype_check(workspace):
    settings, tmp = workspace
    phylo = tmp / "results/kenya/phylo"
    (phylo / "A1/clock_sensitivity.json").write_text(
        json.dumps({
            "estimated": {"clock_rate": 0.00054, "root_date": 1970.13},
            "fixed": {"0.0008": {"root_date": 1979.38, "root_shift_years": 9.25}},
        })
    )  # fmt: skip
    (phylo / "subtype_tree").mkdir()
    (phylo / "subtype_tree/subtype_tree_summary.json").write_text(
        json.dumps({"n_queries": 60, "n_references": 55, "min_support": 95.0,
                    "concordance": {"agree": 5, "unresolved": 50, "not_comparable": 5}})
    )  # fmt: skip
    run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert "Clock-rate sensitivity" in html and "+9.2" in html and "1979.38" in html
    assert "Second subtype call" in html and "Unresolved is common" in html
    assert "unresolved 50" in html


def test_limitations_name_the_scope_of_the_build(tmp_path):
    import shutil

    shutil.copytree("docs", tmp_path / "docs")
    africa = limitations(tmp_path, "africa")
    assert any("the Africa epidemic" in i for i in africa)
    assert not any("Kenya" in i and "epidemic" in i for i in africa)
    assert any("the Kenyan epidemic" in i for i in limitations(tmp_path, "kenya"))


def test_early_sample_caveat_is_in_every_report_and_the_example_only_in_africa(tmp_path):
    import shutil

    shutil.copytree("docs", tmp_path / "docs")
    kenya = [i for i in limitations(tmp_path, "kenya") if "before antiretroviral" in i]
    africa = [i for i in limitations(tmp_path, "africa") if "before antiretroviral" in i]
    assert len(kenya) == len(africa) == 1
    assert "1986" not in kenya[0] and "M41L" in africa[0]


def test_underflowed_p_values_are_shown_as_a_bound_not_zero():
    from virotrace.export.report import _p_text

    assert _p_text("0") == "< 1e-300" and _p_text("8.07e-182") == "8.07e-182"


def test_a_root_older_than_hiv_is_marked_not_plausible():
    from virotrace.export.report import _root_text

    assert "not plausible" in _root_text(1711.39) and _root_text(1962.21) == "1962.21"


def test_the_pooled_row_names_the_country_that_dominates_it():
    from virotrace.export.report import _pooled_note

    rows = [
        {"stratum_kind": "country", "drug_class": "NRTI", "stratum": "A", "denominator": "30"},
        {"stratum_kind": "country", "drug_class": "NRTI", "stratum": "B", "denominator": "70"},
    ]
    note = _pooled_note(rows, ["NRTI"])
    assert "B alone is 70%" in note and _pooled_note(rows[:1], ["NRTI"]) == ""


def test_the_report_ends_with_references_before_provenance(tmp_path):
    from virotrace.export.report import REFERENCES

    assert any("PLoS ONE" in r for r in REFERENCES) and len(REFERENCES) >= 8
