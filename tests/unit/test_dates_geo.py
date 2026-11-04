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
