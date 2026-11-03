"""Configuration contract: valid configs load, known-bad configs are rejected."""

from __future__ import annotations

import pytest

from virotrace.errors import ConfigError
from virotrace.settings import deep_merge, load_settings


def test_base_config_is_valid(config_path):
    settings = load_settings(config_path)
    assert settings.scope == "kenya"
    assert settings.schema_version == 1
    assert len(settings.config_hash()) == 12


@pytest.mark.parametrize("profile", ["dev", "kenya", "africa"])
def test_every_profile_is_valid(config_path, profile):
    settings = load_settings(config_path, profile=profile)
    assert settings.scope in settings.scopes


def test_dev_profile_lowers_suppression_and_warns(config_path):
    settings = load_settings(config_path, profile="dev")
    assert settings.report.min_stratum_n == 1
    assert any("min_stratum_n" in warning for warning in settings.warnings())


def test_profile_overrides_merge_without_restating_siblings(config_path):
    base = load_settings(config_path)
    dev = load_settings(config_path, profile="dev")
    # dev.yaml sets only max_per_stratum under subsample.
    assert dev.subsample.max_per_stratum == 5
    assert dev.subsample.year_bin_width == base.subsample.year_bin_width


def test_deep_merge_does_not_mutate_inputs():
    base = {"a": {"b": 1, "c": 2}}
    override = {"a": {"b": 9}}
    merged = deep_merge(base, override)
    assert merged == {"a": {"b": 9, "c": 2}}
    assert base == {"a": {"b": 1, "c": 2}}


# ------------------------------------------------------------------ known-bad configs
def test_subsample_bounds_inverted_is_rejected(config_path):
    with pytest.raises(ConfigError, match="min_per_stratum"):
        load_settings(config_path, overrides={"subsample": {"min_per_stratum": 99}})


def test_undefined_scope_is_rejected(config_path):
    with pytest.raises(ConfigError, match="not defined in scopes"):
        load_settings(config_path, overrides={"scope": "mars"})


def test_limitations_cannot_be_switched_off(config_path):
    """Not configurable by design -- docs/06-configuration.md section 6."""
    with pytest.raises(ConfigError, match="include_limitations"):
        load_settings(config_path, overrides={"report": {"include_limitations": False}})


def test_coverage_requirement_consistency_is_enforced(config_path):
    with pytest.raises(ConfigError, match="min_gene_coverage"):
        load_settings(
            config_path,
            overrides={"qc": {"min_gene_coverage": 0.0}, "drm": {"require_gene_coverage": True}},
        )


def test_non_codon_gene_range_is_rejected(config_path):
    with pytest.raises(ConfigError, match="multiple of 3"):
        load_settings(config_path, overrides={"reference": {"genes": {"rt": [2550, 4230]}}})


def test_forced_admin1_warns_loudly(config_path):
    settings = load_settings(
        config_path, overrides={"phylo": {"traits": {"admin1": {"enabled": True}}}}
    )
    assert any("overstates" in warning for warning in settings.warnings())


def test_disabled_hypermutation_warns(config_path):
    settings = load_settings(config_path, overrides={"qc": {"hypermutation": {"enabled": False}}})
    assert any("APOBEC" in warning for warning in settings.warnings())


# ------------------------------------------------------------------ determinism
def test_stage_seeds_are_stable_and_distinct(config_path):
    settings = load_settings(config_path)
    assert settings.stage_seed("tree") == settings.stage_seed("tree")
    assert settings.stage_seed("tree") != settings.stage_seed("subsample")


def test_config_hash_changes_with_a_parameter(config_path):
    base = load_settings(config_path)
    changed = load_settings(config_path, overrides={"subsample": {"max_per_stratum": 41}})
    assert base.config_hash() != changed.config_hash()
