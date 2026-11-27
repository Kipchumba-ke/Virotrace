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
            "sdrm_cluster_odds_ratio_adjusted_for_study":
                {"odds_ratio": 1.12, "ci_low": 0.85, "ci_high": 1.48,
                 "n_studies": 143, "n_studies_informative": 52},
            "sensitivity_excluding_studies": {
                "excluded_studies": ["S9"], "n_sequences_removed": 7,
                "sdrm_cluster_odds_ratio_adjusted_for_study":
                    {"odds_ratio": 1.22, "ci_low": 0.92, "ci_high": 1.63,
                     "n_studies": 141, "n_studies_informative": 51}},
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
    assert "1.12 (95% interval 0.85 to 1.48)" in html and "52 of 143 studies" in html
    assert "S9 removed" in html and "1.22 (95% interval 0.92 to 1.63)" in html
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
    assert "no bootstrap support values" in html and "95% bootstrap" not in html
    summary = phylo / "subtype_tree/subtype_tree_summary.json"
    summary.write_text(
        json.dumps({"n_queries": 60, "n_references": 55, "min_support": 95.0,
                    "support_used": True,
                    "concordance": {"agree": 5, "unresolved": 50, "not_comparable": 5}})
    )  # fmt: skip
    run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert "95% bootstrap support" in html and "no bootstrap support values" not in html


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


def test_tree_set_section_names_tips_dropped_by_the_neighbour_screen(workspace):
    settings, tmp = workspace
    path = tmp / "results/kenya/qc/subsample_summary.json"
    path.write_text(
        json.dumps({"waterfall": {"retained": 60, "after_neighbour_screen": 58},
                    "year_span": [2008, 2020], "by_subtype": {"A1": 58}, "by_year": {"2010": 58},
                    "neighbour_screen": {"threshold": 0.86, "removed": [
                        {"strain": "KP223764.1", "subtype": "B", "identity": 0.835},
                        {"strain": "KT183251.1", "subtype": "B", "identity": 0.846}]}})
    )  # fmt: skip
    run_report(settings, repo_root=tmp)
    html = (tmp / "reports/kenya/index.html").read_text()
    assert "58 individuals with strictly dated sequences" in html
    assert "2 sequences were left out of the tree" in html and "KP223764.1, KT183251.1" in html
    assert "86%" in html and "remain in the resistance prevalence" in html


def test_site_is_one_page_with_aggregate_data_palette_and_every_limitation(workspace):
    from virotrace.export.site import run_site

    settings, tmp = workspace
    (tmp / "auspice").mkdir()
    (tmp / "auspice/virotrace_kenya.json").write_text("{}")
    (tmp / "auspice/virotrace_kenya_A1.json").write_text("{}")
    summary = run_site(settings, repo_root=tmp, base_url="https://example.org/virotrace/")
    page = (tmp / "reports/kenya/site/index.html").read_text()
    data = json.loads((tmp / "reports/kenya/site/data.json").read_text())
    assert summary["n_trees"] == 2 and summary["site"].endswith("site/index.html")
    assert (tmp / "reports/kenya/site/trees/virotrace_kenya_A1.json").exists()
    assert {t["label"] for t in data["trees"]} == {"pooled", "A1"}
    assert data["trees"][0]["nextstrain"].startswith("https://nextstrain.org/fetch/example.org/")
    # the page carries the palette, the headline numbers and the whole limitations list
    assert "#101D2D" in page and "#13B8A6" in page and "#8B7BE8" in page
    # the in-page tree viewer reads the same Auspice JSON the trees list points at
    assert "function treeViewer()" in page and "Tree viewer" in page and "fetch(t.file)" in page
    assert data["trees"][0]["file"].startswith("trees/")
    assert len(data["limitations"]) == len(limitations(tmp)) and "Not a clinical tool" in page
    assert any(r["proportion"] == 0.25 and r["denominator"] == 100 for r in data["prevalence"])
    assert data["clusters"]["sdrm_cluster_odds_ratio_adjusted_for_study"]["odds_ratio"] == 1.12
    # one data block and one script: nothing in the data can end the block early
    assert (
        page.count("</script>") == 2
        and "http://" not in page.split('<script id="virotrace-data"')[0]
    )
    # aggregate only: no per-sequence identifiers leave the pipeline
    assert "strain" not in data["prevalence"][0] and "accession" not in json.dumps(data)


def test_site_escapes_script_terminators_inside_the_data(workspace):
    from virotrace.export.site import collect, render

    settings, tmp = workspace
    data = collect(settings, tmp)
    data["limitations"].append("bad </script><script>alert(1)</script> <!-- text")
    page = render(data)
    assert page.count("</script>") == 2 and "<script>alert(1)" not in page


def test_site_without_prevalence_names_the_stage(workspace):
    from virotrace.export.site import run_site

    settings, tmp = workspace
    (tmp / "results/kenya/drm/prevalence.tsv").unlink()
    with pytest.raises(ContractError, match="drm prevalence"):
        run_site(settings, repo_root=tmp)


def test_country_names_from_genbank_map_to_the_outline_and_unknowns_are_reported():
    from virotrace.export.site import country_ids

    names = ["Kenya", "Cote d'Ivoire", "Democratic Republic of the Congo", "Congo",
             "Tanzania", "Swaziland", "Gambia", "Atlantis"]  # fmt: skip
    found, missing = country_ids(names)
    assert found["Kenya"] == "KEN" and found["Cote d'Ivoire"] == "CIV"
    assert found["Democratic Republic of the Congo"] == "COD" and found["Congo"] == "COG"
    assert found["Tanzania"] == "TZA" and found["Swaziland"] == "SWZ" and found["Gambia"] == "GMB"
    assert missing == ["Atlantis"]


def test_site_reports_annotation_coverage_by_context_and_embeds_a_map(workspace):
    from virotrace.export.site import collect

    settings, tmp = workspace
    interim = tmp / "interim/kenya"
    write_tsv(interim / "qc.passed.tsv", [{"strain": s, "genes": "pol", "hypermutated": ""}
                                         for s in "abcde"], ["strain", "genes", "hypermutated"])  # fmt: skip
    write_tsv(
        interim / "metadata.valid.tsv",
        [{"strain": s, "study_key": k} for s, k in zip("abcdef", ["S1", "S1", "S2", "S3", "S4", "S1"], strict=True)],
        ["strain", "study_key"],
    )  # fmt: skip
    (tmp / "study_annotations.tsv").write_text(
        "study_key\ttreatment_context\tnote\nS1\tnaive\t\nS2\tmixed\t\nS3\tunknown\t\n"
    )
    data = collect(settings, tmp)
    an = data["annotation"]
    # f is not QC-passed, so it is not counted; S4 has no annotation
    assert an["n_sequences"] == 5 and an["n_annotated"] == 4
    assert an["by_context"] == {"naive": 2, "mixed": 1, "unknown": 1, "unannotated": 1}
    assert an["n_studies"] == 4 and an["n_studies_annotated"] == 3
    assert data["map"] is None  # the fixture has no country rows
    write_tsv(
        tmp / "results/kenya/drm/prevalence.tsv",
        [_row("country", "Kenya", "NRTI", 25, 100), _row("country", "Atlantis", "NRTI", 1, 50)],
        PREV,
    )
    data = collect(settings, tmp)
    assert data["country_ids"] == {"Kenya": "KEN"} and data["unmapped_countries"] == ["Atlantis"]
    assert {f["id"] for f in data["map"]["features"]} >= {"KEN", "UGA", "ZAF", "NGA"}


def test_site_annotation_coverage_is_absent_without_interim_tables(workspace):
    from virotrace.export.site import collect

    settings, tmp = workspace
    assert collect(settings, tmp)["annotation"] is None
