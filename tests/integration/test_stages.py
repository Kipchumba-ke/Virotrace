"""The acquire and normalise stages end to end on fixtures, plus config resolution."""

from __future__ import annotations

import json

import pytest

from virotrace.acquire.genbank import GenBankAdapter
from virotrace.acquire.snapshot import resolve_snapshot
from virotrace.acquire.stage import build_adapters, run_acquire
from virotrace.curate.geo import countries_for_scope
from virotrace.curate.stage import run_normalise
from virotrace.errors import ConfigError, ContractError
from virotrace.io.tables import read_tsv
from virotrace.settings import load_settings


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("VIROTRACE_PROFILE", raising=False)
    monkeypatch.delenv("VIROTRACE_SCOPE", raising=False)


def _settings(config_path, tmp_path, records_toy_tsv, **overrides):
    base = {
        "paths": {
            "raw": str(tmp_path / "raw"),
            "interim": str(tmp_path / "interim"),
            "processed": str(tmp_path / "processed"),
            "reference": "config/reference",
            "results": str(tmp_path / "results"),
            "auspice": str(tmp_path / "auspice"),
            "reports": str(tmp_path / "reports"),
        },
        "acquire": {"sources": ["fixture"], "fixture_path": str(records_toy_tsv)},
    }
    base.update(overrides)
    return load_settings(config_path, "dev", base)


# ------------------------------------------------------------------ kenya scope
def test_kenya_scope_end_to_end(config_path, tmp_path, records_toy_tsv):
    settings = _settings(config_path, tmp_path, records_toy_tsv)
    snapshot = run_acquire(settings)
    summary = run_normalise(settings)

    assert summary["snapshot_id"] == snapshot.name
    assert summary["n_snapshot"] == 9
    assert summary["n_quarantined"] == 1  # Atlantis
    assert summary["n_out_of_scope"] == 3  # Uganda, Tanzania, South Africa
    assert summary["n_retained"] == 5
    assert summary["n_collapsed_by_dedup"] == 0
    assert summary["corpus"]["n_ken"] == 5

    meta = read_tsv(tmp_path / "interim" / "kenya" / "metadata.valid.tsv")
    assert {r["country_code"] for r in meta} == {"KEN"}
    assert len(meta) == 5

    quarantine = read_tsv(tmp_path / "results" / "kenya" / "qc" / "quarantine.tsv")
    assert [q["accession"] for q in quarantine] == ["XX000001"]

    counts = json.loads((tmp_path / "results" / "kenya" / "qc" / "corpus_counts.json").read_text())
    assert counts["n_retained"] == 5


def test_every_record_is_accounted_for(config_path, tmp_path, records_toy_tsv):
    """Conservation across the whole stage: nothing vanishes between snapshot and outputs."""
    settings = _settings(config_path, tmp_path, records_toy_tsv)
    run_acquire(settings)
    s = run_normalise(settings)
    assert s["n_snapshot"] == (
        s["n_collapsed_by_dedup"] + s["n_quarantined"] + s["n_out_of_scope"] + s["n_retained"]
    )


def test_africa_scope_keeps_regional_context(config_path, tmp_path, records_toy_tsv):
    settings = _settings(config_path, tmp_path, records_toy_tsv, scope="africa")
    run_acquire(settings)
    summary = run_normalise(settings)
    assert summary["n_retained"] == 8  # everything African; Atlantis still quarantined
    assert summary["n_out_of_scope"] == 0
    assert (tmp_path / "interim" / "africa" / "metadata.valid.tsv").exists()


# --------------------------------------------------------------------- pointers
def test_acquire_writes_a_scope_pointer_that_normalise_follows(
    config_path, tmp_path, records_toy_tsv
):
    settings = _settings(config_path, tmp_path, records_toy_tsv)
    snapshot = run_acquire(settings)
    pointer = tmp_path / "raw" / "kenya.current"
    assert pointer.read_text().strip() == snapshot.name
    assert resolve_snapshot(tmp_path / "raw", "kenya") == snapshot


def test_pinned_snapshot_id_wins_over_the_pointer(config_path, tmp_path, records_toy_tsv):
    settings = _settings(config_path, tmp_path, records_toy_tsv)
    run_acquire(settings)
    with pytest.raises(ContractError, match="pinned snapshot"):
        resolve_snapshot(tmp_path / "raw", "kenya", pinned_id="2020-01-01_deadbe")


def test_stale_pointer_is_a_loud_error(tmp_path):
    (tmp_path / "kenya.current").write_text("2020-01-01_gone\n")
    with pytest.raises(ContractError, match="no manifest"):
        resolve_snapshot(tmp_path, "kenya")


# ------------------------------------------------------------------ adapter wiring
def test_genbank_source_builds_a_kenya_query(config_path, tmp_path, records_toy_tsv):
    settings = _settings(config_path, tmp_path, records_toy_tsv, acquire={"sources": ["genbank"]})
    (adapter,) = build_adapters(settings)
    assert isinstance(adapter, GenBankAdapter)
    assert adapter.term == "txid11676[Organism] AND Kenya[Country] AND 1000:10000[SLEN]"


def test_africa_scope_builds_an_or_group_of_countries(config_path, tmp_path, records_toy_tsv):
    settings = _settings(
        config_path, tmp_path, records_toy_tsv, scope="africa", acquire={"sources": ["genbank"]}
    )
    (adapter,) = build_adapters(settings)
    assert "(Algeria[Country] OR Angola[Country]" in adapter.term
    assert "Kenya[Country]" in adapter.term


def test_capped_trial_run_has_its_own_snapshot_identity(config_path, tmp_path, records_toy_tsv):
    settings = _settings(config_path, tmp_path, records_toy_tsv, acquire={"sources": ["genbank"]})
    (capped,) = build_adapters(settings, max_records=20)
    (full,) = build_adapters(settings)
    assert capped.query_description() != full.query_description()
    assert capped.query_description().endswith("max_records=20")


def test_unimplemented_and_unknown_sources_fail_clearly(config_path, tmp_path, records_toy_tsv):
    for source, message in (("lanl", "not implemented"), ("nope", "unknown acquire source")):
        settings = _settings(config_path, tmp_path, records_toy_tsv, acquire={"sources": [source]})
        with pytest.raises(ConfigError, match=message):
            build_adapters(settings)


def test_countries_for_scope():
    assert countries_for_scope(["KEN"], []) == ["Kenya"]
    east = countries_for_scope([], ["southern_africa"])
    assert east == ["Botswana", "Eswatini", "Lesotho", "Namibia", "South Africa"]
    assert countries_for_scope(["XXX"], []) == []


# ------------------------------------------------------ config resolution via env
def test_env_variables_resolve_like_cli_arguments(config_path, monkeypatch):
    monkeypatch.setenv("VIROTRACE_PROFILE", "dev")
    monkeypatch.setenv("VIROTRACE_SCOPE", "africa")
    settings = load_settings(config_path)
    assert settings.scope == "africa"  # env scope applied after the profile's scope
    assert settings.report.min_stratum_n == 1  # dev profile applied


def test_explicit_arguments_beat_the_environment(config_path, monkeypatch):
    monkeypatch.setenv("VIROTRACE_PROFILE", "dev")
    settings = load_settings(config_path, profile="africa")
    assert settings.report.min_stratum_n == 15  # africa profile, not dev


# ------------------------------------------------------------------ clone policy
def _fixture_with_clones(records_toy_tsv, tmp_path):
    """A copy of the toy fixture with a /clone tag on two Kenyan records."""
    lines = records_toy_tsv.read_text().splitlines()
    header, body = lines[0], lines[1:]
    out = [header + "\tclone"]
    for line in body:
        accession = line.split("\t")[0]
        out.append(f"{line}\t{'c1' if accession in ('KE000001', 'KE000002') else ''}")
    path = tmp_path / "records_clones.tsv"
    path.write_text("\n".join(out) + "\n")
    return path


def test_clone_policy_excludes_from_downstream_but_keeps_the_rows(
    config_path, tmp_path, records_toy_tsv
):
    fixture = _fixture_with_clones(records_toy_tsv, tmp_path)
    settings = _settings(config_path, tmp_path, fixture)
    run_acquire(settings)
    summary = run_normalise(settings)

    assert summary["clone_policy"] == "exclude"
    assert summary["n_clone_derived"] == 2
    assert summary["n_retained"] == 5  # conservation is untouched by the policy
    assert summary["n_after_clone_policy"] == 3
    assert summary["corpus_after_clone_policy"]["n_total"] == 3

    meta = read_tsv(tmp_path / "interim" / "kenya" / "metadata.valid.tsv")
    assert len(meta) == 5  # clone-derived rows stay in the table for the DRM layer
    flagged = {r["strain"] for r in meta if r["clone_policy_excluded"] == "true"}
    assert flagged == {"KE000001", "KE000002"}


def test_keep_all_policy_flags_nothing_for_exclusion(config_path, tmp_path, records_toy_tsv):
    fixture = _fixture_with_clones(records_toy_tsv, tmp_path)
    settings = _settings(config_path, tmp_path, fixture, qc={"clone_policy": "keep_all"})
    run_acquire(settings)
    summary = run_normalise(settings)

    assert summary["n_clone_derived"] == 2  # still reported
    assert summary["n_after_clone_policy"] == 5
    meta = read_tsv(tmp_path / "interim" / "kenya" / "metadata.valid.tsv")
    assert all(r["clone_policy_excluded"] == "false" for r in meta)
