"""`phylo traits`: command, admin1 gating and outputs, with augur injected."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from virotrace.errors import ContractError
from virotrace.io.tables import write_tsv
from virotrace.phylo.traits import resolved_fraction, run_traits
from virotrace.settings import load_settings


def _settings(config_path, tmp_path, admin1="auto", min_frac=0.3):
    return load_settings(
        config_path,
        "dev",
        {
            "paths": {
                "raw": str(tmp_path / "raw"), "interim": str(tmp_path / "interim"),
                "processed": str(tmp_path / "processed"), "reference": str(tmp_path),
                "results": str(tmp_path / "results"), "auspice": str(tmp_path / "auspice"),
                "reports": str(tmp_path / "reports"),
            },
            "phylo": {
                "traits": {
                    "enabled": True,
                    "admin1": {"enabled": admin1, "min_resolved_frac": min_frac},
                }
            },
        },
    )  # fmt: skip


def _inputs(tmp_path, resolved):
    phylo = tmp_path / "results/kenya/phylo"
    phylo.mkdir(parents=True)
    (phylo / "tree.refined.nwk").write_text("((A:1,B:1):1,C:2);\n")
    rows = [
        {"strain": s, "country_code": "KE", "region": "East Africa" if s != "C" else "",
         "admin1": "Nairobi", "admin1_resolved": r}
        for s, r in zip("ABC", resolved, strict=True)
    ]  # fmt: skip
    write_tsv(tmp_path / "processed/kenya/metadata.final.tsv", rows)
    return phylo


def _fake_augur(command):
    out = command[command.index("--output-node-data") + 1]
    cols = command[command.index("--columns") + 1 : command.index("--confidence")]
    Path(out).write_text(json.dumps({"models": {c: {} for c in cols}, "nodes": {}}))


def test_traits_runs_augur_on_the_refined_tree_and_marks_missing_values(config_path, tmp_path):
    settings = _settings(config_path, tmp_path)
    phylo = _inputs(tmp_path, ["False", "False", "False"])
    calls = []

    def runner(command):
        calls.append(command)
        _fake_augur(command)

    out = run_traits(settings, repo_root=".", runner=runner)
    cmd = calls[0]
    assert cmd[:2] == ["augur", "traits"] and "--confidence" in cmd
    assert out["columns"] == ["country_code", "region"]  # admin1 gated out: 0% resolved
    assert out["admin1"]["reconstructed"] is False and out["n_missing"]["region"] == 1
    assert "C\tKE\t?" in (phylo / "metadata.traits.tsv").read_text()
    assert json.loads((phylo / "traits_summary.json").read_text())["n_tips"] == 3


def test_admin1_is_reconstructed_only_when_enough_is_resolved(config_path, tmp_path):
    _inputs(tmp_path, ["True", "True", "False"])
    got = run_traits(_settings(config_path, tmp_path), repo_root=".", runner=_fake_augur)
    assert got["columns"][-1] == "admin1" and got["admin1"]["resolved_fraction"] == 0.6667
    off = run_traits(_settings(config_path, tmp_path, False), repo_root=".", runner=_fake_augur)
    assert "admin1" not in off["columns"]
    strict = run_traits(
        _settings(config_path, tmp_path, "auto", 0.9), repo_root=".", runner=_fake_augur
    )
    assert "admin1" not in strict["columns"]


def test_resolved_fraction_handles_empty_and_text_flags():
    assert resolved_fraction([]) == 0.0
    assert resolved_fraction([{"admin1_resolved": "True"}, {"admin1_resolved": "False"}]) == 0.5


def test_missing_inputs_and_no_output_are_named(config_path, tmp_path):
    settings = _settings(config_path, tmp_path)
    with pytest.raises(ContractError, match="phylo refine"):
        run_traits(settings, repo_root=".", runner=_fake_augur)
    _inputs(tmp_path, ["False"] * 3)
    with pytest.raises(ContractError, match="no node data"):
        run_traits(settings, repo_root=".", runner=lambda c: None)
