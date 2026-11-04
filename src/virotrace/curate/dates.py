"""Collection-date parsing with Nextstrain-style ambiguity notation.

Precision is never invented: a year-only date stays ``YYYY-XX-XX`` and is marked
``year`` precision, which decides whether the record can enter the time tree.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

MONTHS = ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

_YEAR = re.compile(r"^(?P<y>\d{4})$")
_ISO_DAY = re.compile(r"^(?P<y>\d{4})[-/](?P<m>\d{1,2})[-/](?P<d>\d{1,2})$")
_ISO_MONTH = re.compile(r"^(?P<y>\d{4})[-/](?P<m>\d{1,2})$")
_GB_DAY = re.compile(r"^(?P<d>\d{1,2})[-\s](?P<m>[A-Za-z]{3})[a-z]*[-\s](?P<y>\d{4})$")
_GB_MONTH = re.compile(r"^(?P<m>[A-Za-z]{3})[a-z]*[-\s](?P<y>\d{4})$")


@dataclass(frozen=True, slots=True)
class ParsedDate:
    date: str  # YYYY-MM-DD | YYYY-MM-XX | YYYY-XX-XX | "" when unknown
    precision: str  # day | month | year | unknown
    year: int | None


UNKNOWN = ParsedDate("", "unknown", None)


def _month(token: str) -> int | None:
    token = token.lower()[:3]
    return MONTHS.index(token) + 1 if token in MONTHS else None


def parse_collection_date(
    raw: str | None, *, min_year: int = 1985, max_year: int = 2026
) -> ParsedDate:
    """Parse a source date string. Out-of-range or malformed values become unknown."""
    text = (raw or "").strip()
    if not text:
        return UNKNOWN

    year: int | None = None
    month: int | None = None
    day: int | None = None

    if m := _YEAR.match(text):
        year = int(m["y"])
    elif m := _ISO_DAY.match(text):
        year, month, day = int(m["y"]), int(m["m"]), int(m["d"])
    elif m := _ISO_MONTH.match(text):
        year, month = int(m["y"]), int(m["m"])
    elif m := _GB_DAY.match(text):
        year, month, day = int(m["y"]), _month(m["m"]), int(m["d"])
    elif m := _GB_MONTH.match(text):
        year, month = int(m["y"]), _month(m["m"])
    else:
        return UNKNOWN

    if year is None or not min_year <= year <= max_year:
        return UNKNOWN
    if month is not None and not 1 <= month <= 12:
        return UNKNOWN
    if day is not None and not 1 <= day <= 31:
        return UNKNOWN

    if month is None:
        return ParsedDate(f"{year:04d}-XX-XX", "year", year)
    if day is None:
        return ParsedDate(f"{year:04d}-{month:02d}-XX", "month", year)
    return ParsedDate(f"{year:04d}-{month:02d}-{day:02d}", "day", year)
