"""Deterministic stratified subsampling.

Random downsampling would bias composition toward whichever country or study
deposited most sequences. Stratifying on (country, year bin, subtype) and ranking
within a stratum means reducing n costs precision without shifting composition.

Determinism: inputs are sorted, priority keys are explicit, and ties break on
``strain``. The same inputs and seed always give the same output -- property-tested.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

DATE_PRECISION_RANK = {"day": 0, "month": 1, "year": 2, "unknown": 3}


@dataclass(frozen=True, slots=True)
class StratumSummary:
    key: tuple[Any, ...]
    n_available: int
    n_retained: int


def year_bin(year: int, width: int) -> int:
    """Left edge of the bin containing ``year``."""
    if width < 1:
        raise ValueError("width must be >= 1")
    return (year // width) * width


def _priority_key(record: Mapping[str, Any], priority: Sequence[str]) -> tuple[Any, ...]:
    parts: list[Any] = []
    for field in priority:
        value = record.get(field)
        if field == "date_precision":
            parts.append(DATE_PRECISION_RANK.get(str(value), 99))
        elif field == "gene_coverage":
            parts.append(-float(value if value is not None else 0.0))  # higher is better
        elif field == "ambiguity_frac":
            parts.append(float(value if value is not None else 1.0))  # lower is better
        else:
            parts.append(str(value))
    parts.append(str(record.get("strain", "")))  # deterministic tie-break
    return tuple(parts)


def stratified_subsample(
    records: Iterable[Mapping[str, Any]],
    *,
    strata: Sequence[str],
    max_per_stratum: int,
    min_per_stratum: int = 1,
    priority: Sequence[str] = ("date_precision", "gene_coverage", "ambiguity_frac"),
    year_bin_width: int = 2,
) -> tuple[list[str], list[StratumSummary]]:
    """Return ``(retained_strains, stratum_summaries)``.

    Idempotent for a given input set; never exceeds ``max_per_stratum`` per stratum;
    never drops a non-empty stratum below ``min_per_stratum``.
    """
    if min_per_stratum > max_per_stratum:
        raise ValueError("min_per_stratum exceeds max_per_stratum")

    buckets: dict[tuple[Any, ...], list[Mapping[str, Any]]] = {}
    for record in records:
        key_parts: list[Any] = []
        for field in strata:
            if field == "year_bin":
                key_parts.append(year_bin(int(record["year"]), year_bin_width))
            else:
                key_parts.append(record.get(field))
        buckets.setdefault(tuple(key_parts), []).append(record)

    retained: list[str] = []
    summaries: list[StratumSummary] = []
    for key in sorted(buckets, key=lambda k: tuple(str(part) for part in k)):
        members = sorted(buckets[key], key=lambda r: _priority_key(r, priority))
        take = max(min(max_per_stratum, len(members)), min(min_per_stratum, len(members)))
        chosen = members[:take]
        retained.extend(str(record["strain"]) for record in chosen)
        summaries.append(StratumSummary(key, len(members), len(chosen)))

    return sorted(retained), summaries
