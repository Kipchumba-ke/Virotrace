"""Stratified subsampling: determinism, caps, and composition preservation."""

from __future__ import annotations

import pytest

from virotrace.curate.subsample import stratified_subsample, year_bin


def test_year_bin_edges():
    assert year_bin(2015, 2) == 2014
    assert year_bin(2016, 2) == 2016
    assert year_bin(2016, 5) == 2015
    with pytest.raises(ValueError):
        year_bin(2016, 0)


def test_cap_is_never_exceeded(toy_records):
    _, summaries = stratified_subsample(
        toy_records,
        strata=["country_code", "year_bin", "subtype_called"],
        max_per_stratum=1,
    )
    assert all(summary.n_retained <= 1 for summary in summaries)


def test_non_empty_stratum_is_never_dropped(toy_records):
    retained, summaries = stratified_subsample(
        toy_records,
        strata=["country_code", "year_bin", "subtype_called"],
        max_per_stratum=1,
        min_per_stratum=1,
    )
    assert all(summary.n_retained >= 1 for summary in summaries)
    # Every stratum present in the input is represented in the output.
    assert len(summaries) == len(
        {(r["country_code"], year_bin(r["year"], 2), r["subtype_called"]) for r in toy_records}
    )
    assert len(retained) == len(summaries)


def test_is_idempotent_and_order_independent(toy_records):
    first, _ = stratified_subsample(toy_records, strata=["country_code"], max_per_stratum=2)
    second, _ = stratified_subsample(toy_records, strata=["country_code"], max_per_stratum=2)
    shuffled, _ = stratified_subsample(
        list(reversed(toy_records)), strata=["country_code"], max_per_stratum=2
    )
    assert first == second == shuffled


def test_priority_prefers_better_dated_and_better_covered_records(toy_records):
    """Reducing n should cost precision, not shift composition or quality."""
    retained, _ = stratified_subsample(
        toy_records,
        strata=["country_code"],
        max_per_stratum=1,
    )
    # Within KEN, KE_001 has day precision and the highest coverage.
    assert "KE_001" in retained
    assert "KE_005" not in retained  # unknown date precision, lowest coverage


def test_rejects_inverted_bounds(toy_records):
    with pytest.raises(ValueError):
        stratified_subsample(
            toy_records, strata=["country_code"], max_per_stratum=1, min_per_stratum=5
        )


def test_minority_subtype_survives_capping():
    """Random downsampling would lose the rare stratum; stratification must not."""
    records = [
        {
            "strain": f"A{i}",
            "country_code": "KEN",
            "year": 2016,
            "subtype_called": "A1",
            "date_precision": "day",
            "gene_coverage": 0.9,
            "ambiguity_frac": 0.0,
        }
        for i in range(200)
    ] + [
        {
            "strain": "D1",
            "country_code": "KEN",
            "year": 2016,
            "subtype_called": "D",
            "date_precision": "day",
            "gene_coverage": 0.9,
            "ambiguity_frac": 0.0,
        }
    ]
    retained, _ = stratified_subsample(
        records, strata=["country_code", "year_bin", "subtype_called"], max_per_stratum=5
    )
    assert "D1" in retained
    assert len(retained) == 6
