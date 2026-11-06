"""Auspice export: tip metadata, colourings and the augur commands, with an injected runner."""

from __future__ import annotations

import json

import pytest

from virotrace.errors import ContractError
from virotrace.export.auspice import auspice_config, run_export, tip_metadata
from virotrace.io.fasta import write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings


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
    processed, results = tmp_path / "processed/kenya", tmp_path / "results/kenya"
    strains = ["s1", "s2", "s3", "s4"]
    write_tsv(
        processed / "metadata.final.tsv",
        [
            {"strain": s, "country": "Kenya", "region": "eastern_africa", "year": "2012",
             "date": "2012-05-03", "date_precision": "day",
             "subtype_called": "A1" if s != "s4" else "D",
             "study_key": "S1" if s != "s4" else "S2"}
            for s in strains
        ],
    )  # fmt: skip
    write_tsv(
        results / "phylo/clusters.tsv",
        [{"strain": s, "cluster_id": "CL001" if s in ("s1", "s2") else "", "cluster_size": 2 if s in ("s1", "s2") else 1}
         for s in strains],
    )  # fmt: skip
    write_tsv(
        results / "drm/drm_table.tsv",
        [
            {"strain": "s1", "gene": "RT", "position": 103, "ref_aa": "K", "sdrm_residues": "N"},
            {"strain": "s1", "gene": "RT", "position": 184, "ref_aa": "M", "sdrm_residues": "V"},
            {"strain": "s3", "gene": "PR", "position": 46, "ref_aa": "M", "sdrm_residues": "I"},
        ],
        ["strain", "gene", "position", "ref_aa", "sdrm_residues"],
    )  # fmt: skip
    (tmp_path / "study_annotations.tsv").write_text("study_key\ttreatment_context\nS1\tnaive\n")
    write_fasta(processed / "aln.tree.fasta", [(s, "ACGT" * 5) for s in strains], wrap=0)
    write_fasta(processed / "aln.tree.A1.fasta", [(s, "ACGT" * 5) for s in strains[:3]], wrap=0)
    for folder, time_tree in ((results / "phylo", False), (results / "phylo/A1", True)):
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "tree.refined.nwk").write_text("(s1:1,s2:1,s3:1);\n")
        (folder / "branch_lengths.json").write_text("{}")
        (folder / "temporal_gate.json").write_text(json.dumps({"time_tree": time_tree}))
    return settings, tmp_path


def test_tip_metadata_has_clusters_treatment_context_and_sdrm_by_class(workspace):
    settings, _ = workspace
    rows = {r["strain"]: r for r in tip_metadata(settings, ".")}
    assert rows["s1"]["cluster"] == "CL001" and rows["s3"]["cluster"] == "unclustered"
    assert rows["s1"]["treatment_context"] == "naive"
    assert rows["s4"]["treatment_context"] == "unannotated"
    assert (rows["s1"]["sdrm_nnrti"], rows["s1"]["sdrm_nrti"], rows["s1"]["sdrm_pi"]) == (
        "yes", "yes", "no",
    )  # fmt: skip
    assert rows["s3"]["sdrm_pi"] == "yes" and rows["s3"]["sdrm_mutations"] == "PR:M46I"
    assert rows["s2"]["sdrm_any"] == "no" and rows["s2"]["sdrm_mutations"] == "none"
    assert rows["s1"]["sdrm_mutations"] == "RT:K103N,RT:M184V"
    assert rows["s1"]["date"] == "2012-05-03"


def test_export_runs_ancestral_then_export_for_the_pooled_and_each_subtype_tree(workspace):
    settings, tmp = workspace
    calls = []
    summary = run_export(settings, repo_root=".", runner=calls.append)
    assert [c[1] for c in calls] == ["ancestral", "export", "ancestral", "export"]
    pooled, a1 = summary["datasets"]["pooled"], summary["datasets"]["A1"]
    assert "divergence tree" in pooled["title"] and not pooled["time_tree"]
    assert "subtype A1, time tree" in a1["title"] and a1["time_tree"]
    assert pooled["dataset"].endswith("virotrace_kenya.json")
    assert a1["dataset"].endswith("virotrace_kenya_A1.json")
    export = calls[3]
    assert export[export.index("--output") + 1] == a1["dataset"]
    assert any(str(a).endswith("branch_lengths.json") for a in export)
    assert str(tmp / "results/kenya/phylo/A1/nt_muts.json") in export
    assert str(tmp / "processed/kenya/aln.tree.A1.fasta") in calls[2]
    config = json.loads((tmp / "results/kenya/phylo/A1/auspice_config.json").read_text())
    assert {c["key"] for c in config["colorings"]} >= {"subtype", "cluster", "study"}
    assert read_tsv(tmp / "results/kenya/phylo/auspice_metadata.tsv")[0]["strain"] == "s1"


def test_colourings_only_name_columns_that_exist():
    config = auspice_config("t", ["strain", "subtype", "study"])
    assert [c["key"] for c in config["colorings"]] == ["subtype", "study"]


def test_missing_inputs_name_the_stage_to_run(workspace):
    settings, tmp = workspace
    (tmp / "results/kenya/phylo/A1/tree.refined.nwk").unlink()
    with pytest.raises(ContractError, match="phylo subtrees"):
        run_export(settings, repo_root=".", runner=lambda c: None)
    (tmp / "results/kenya/phylo/clusters.tsv").unlink()
    with pytest.raises(ContractError, match="phylo clusters"):
        run_export(settings, repo_root=".", runner=lambda c: None)
