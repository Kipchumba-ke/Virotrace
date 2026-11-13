"""The P1 gate: how many Kenyan sequences can actually support a time tree?

Run this before committing to the Phase 3 phylogenomics build. If dated Kenyan
sequences are too few, that is a finding to report, and the fallback (pool East
Africa for the tree, keep Kenya as a trait) is a different build from the spec's
default.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from typing import Any


def count_corpus(rows: Iterable[Mapping[str, Any]], *, country_code: str = "KEN") -> dict[str, Any]:
    rows = list(rows)
    target = [r for r in rows if r["country_code"] == country_code]
    precision = Counter(r["date_precision"] for r in target)
    everywhere = Counter(r["date_precision"] for r in rows)
    years = sorted({r["year"] for r in target if r["year"] is not None})
    return {
        "n_total": len(rows),
        f"n_{country_code.lower()}": len(target),
        "by_precision": {p: precision.get(p, 0) for p in ("day", "month", "year", "unknown")},
        "n_time_tree_eligible_strict": precision.get("day", 0) + precision.get("month", 0),
        "n_time_tree_eligible_with_year_only": (
            precision.get("day", 0) + precision.get("month", 0) + precision.get("year", 0)
        ),
        # the same two counts over the whole corpus, for a scope wider than the target country
        "all_n_time_tree_eligible_strict": everywhere.get("day", 0) + everywhere.get("month", 0),
        "all_n_time_tree_eligible_with_year_only": (
            everywhere.get("day", 0) + everywhere.get("month", 0) + everywhere.get("year", 0)
        ),
        # Where the tree's dated tips actually sit: a clock needs spread, not just a count.
        "strict_dated_by_year": {
            str(y): n
            for y, n in sorted(
                Counter(
                    r["year"]
                    for r in target
                    if r["date_precision"] in ("day", "month") and r["year"] is not None
                ).items()
            )
        },
        "year_range": (years[0], years[-1]) if years else None,
        "n_distinct_years": len(years),
        "by_region": dict(Counter(r["region"] for r in rows)),
    }
