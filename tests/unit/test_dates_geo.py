"""Date parsing and country normalisation."""

from __future__ import annotations

import pytest

from virotrace.curate.dates import parse_collection_date
from virotrace.curate.geo import normalise_country


@pytest.mark.parametrize(
    ("raw", "date", "precision", "year"),
    [
        ("2015-03-12", "2015-03-12", "day", 2015),
        ("12-Mar-2015", "2015-03-12", "day", 2015),
        ("2016-05", "2016-05-XX", "month", 2016),
        ("Jun-2017", "2017-06-XX", "month", 2017),
        ("2014", "2014-XX-XX", "year", 2014),
        ("", "", "unknown", None),
        (None, "", "unknown", None),
        ("sometime in spring", "", "unknown", None),
        ("2015-13-01", "", "unknown", None),
        ("1850", "", "unknown", None),
        ("2099", "", "unknown", None),
    ],
)
def test_date_parsing(raw, date, precision, year):
    parsed = parse_collection_date(raw)
    assert (parsed.date, parsed.precision, parsed.year) == (date, precision, year)


def test_precision_is_never_invented():
    """A year-only date must stay year precision, not become 1 January."""
    parsed = parse_collection_date("2014")
    assert parsed.precision == "year"
    assert "01-01" not in parsed.date


def test_country_with_admin_hint_keeps_it_unresolved():
    geo = normalise_country("Kenya: Nairobi")
    assert geo is not None
    assert (geo.country_code, geo.region) == ("KEN", "eastern_africa")
    assert geo.admin1_raw == "Nairobi"
    assert geo.admin1_resolved is False  # a hint is never promoted to a resolved county


@pytest.mark.parametrize(
    ("raw", "code", "region"),
    [
        ("Uganda", "UGA", "eastern_africa"),
        ("United Republic of Tanzania", "TZA", "eastern_africa"),
        ("South Africa", "ZAF", "southern_africa"),
        ("Ivory Coast", "CIV", "western_africa"),
        ("Cameroon", "CMR", "middle_africa"),
        ("USA", "USA", "non_africa"),
    ],
)
def test_country_and_region_mapping(raw, code, region):
    geo = normalise_country(raw)
    assert geo is not None
    assert (geo.country_code, geo.region) == (code, region)


def test_unmapped_country_is_none_not_a_guess():
    assert normalise_country("Atlantis") is None
    assert normalise_country("") is None
    assert normalise_country(None) is None


@pytest.mark.parametrize("raw", ["Finland", "Viet Nam: Hanoi", "Japan", "Russian Federation"])
def test_real_countries_outside_africa_are_recognised_not_unmapped(raw):
    # An out-of-scope country is a valid record to count, never a data error.
    result = normalise_country(raw)
    assert result is not None
    assert result.region == "non_africa"


def test_corpus_reports_strict_dated_tips_by_year():
    from virotrace.curate.corpus import count_corpus

    rows = [
        {"country_code": "KEN", "date_precision": p, "year": y, "region": "eastern_africa"}
        for p, y in [
            ("day", 2010),
            ("month", 2010),
            ("month", 2012),
            ("year", 2012),
            ("unknown", None),
        ]
    ]
    assert count_corpus(rows)["strict_dated_by_year"] == {"2010": 2, "2012": 1}


def test_corpus_counts_dated_rows_over_the_whole_corpus_as_well_as_the_target_country():
    from virotrace.curate.corpus import count_corpus

    rows = [
        {"country_code": c, "date_precision": p, "year": 2010, "region": "x"}
        for c, p in (("KEN", "day"), ("UGA", "day"), ("UGA", "month"), ("UGA", "year"))
    ]
    out = count_corpus(rows)
    assert out["n_time_tree_eligible_strict"] == 1
    assert out["all_n_time_tree_eligible_strict"] == 3
    assert out["all_n_time_tree_eligible_with_year_only"] == 4
