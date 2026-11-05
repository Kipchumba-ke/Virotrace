"""DRM prevalence estimation with per-class denominators and Wilson intervals.

Two decisions here carry the credibility of the headline numbers:

1. **Per-class denominators.** A sequence with no integrase coverage contributes to
   no INSTI denominator. Using one global n would turn missing data into apparent
   susceptibility.
2. **Wilson intervals.** Wald intervals misbehave at small n and near 0%, which is
   exactly where surveillance estimates live. Wilson is not configurable away.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from typing import Protocol

SUPPRESSED = "n_too_small"


class ClassCoverage(Protocol):
    """What this module needs from a sequence's calls.

    Structural typing keeps prevalence free of an import from ``drm.mutations``,
    which would otherwise be circular.
    """

    def covers_class(self, drug_class: str, min_fraction: float) -> bool: ...


@dataclass(frozen=True, slots=True)
class PrevalenceRow:
    stratum: tuple[str, ...]
    drug_class: str
    numerator: int
    denominator: int
    proportion: float | None
    ci_low: float | None
    ci_high: float | None
    status: str  # "estimated" | "n_too_small"

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _z_for(level: float) -> float:
    """Two-sided normal quantile. Exact for the levels that matter, else Acklam-free."""
    table = {0.90: 1.6448536269514722, 0.95: 1.959963984540054, 0.99: 2.5758293035489004}
    if level in table:
        return table[level]
    # Inverse error function via Newton iteration -- adequate and dependency-free.
    target = level
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        coverage = math.erf(mid / math.sqrt(2))
        if coverage < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def wilson_ci(successes: int, total: int, level: float = 0.95) -> tuple[float, float]:
    """Wilson score interval. Bounds are clamped to [0, 1] and contain the estimate."""
    if total <= 0:
        raise ValueError("total must be positive")
    if not 0 <= successes <= total:
        raise ValueError(f"successes ({successes}) outside 0..{total}")

    z = _z_for(level)
    p = successes / total
    denom = 1 + z**2 / total
    centre = (p + z**2 / (2 * total)) / denom
    margin = (z / denom) * math.sqrt(p * (1 - p) / total + z**2 / (4 * total**2))

    # At p=0 the lower bound is analytically exactly 0, and at p=1 the upper bound
    # is analytically exactly 1. In floating point both come out a few ulps off
    # (1e-19 rather than 0.0), which breaks the "interval contains the estimate"
    # invariant at the boundaries. Snap them rather than clamp.
    low = 0.0 if successes == 0 else max(0.0, centre - margin)
    high = 1.0 if successes == total else min(1.0, centre + margin)
    return low, high


def estimate_prevalence(
    *,
    numerator: int,
    denominator: int,
    drug_class: str,
    stratum: tuple[str, ...] = (),
    min_stratum_n: int = 10,
    ci_level: float = 0.95,
) -> PrevalenceRow:
    """One stratum's prevalence, suppressed if the denominator is too small.

    Suppression is reported as a row with ``status='n_too_small'`` rather than
    dropped, so the report can show that the stratum existed and was not estimable.
    An unsuppressed 100% from n=2 is the classic surveillance-dashboard error.
    """
    if numerator > denominator:
        raise ValueError(
            f"numerator ({numerator}) exceeds denominator ({denominator}) for {drug_class}"
        )
    if denominator < min_stratum_n:
        return PrevalenceRow(
            stratum=stratum,
            drug_class=drug_class,
            numerator=numerator,
            denominator=denominator,
            proportion=None,
            ci_low=None,
            ci_high=None,
            status=SUPPRESSED,
        )
    low, high = wilson_ci(numerator, denominator, ci_level)
    return PrevalenceRow(
        stratum=stratum,
        drug_class=drug_class,
        numerator=numerator,
        denominator=denominator,
        proportion=numerator / denominator,
        ci_low=low,
        ci_high=high,
        status="estimated",
    )


def cochran_armitage_trend(
    numerators: Sequence[int], denominators: Sequence[int], scores: Sequence[float] | None = None
) -> tuple[float, float]:
    """Cochran-Armitage test for trend across ordered strata.

    Returns ``(z, p_two_sided)``. Interpretation caveat, carried into the report:
    sampling composition changes over time, so a significant trend may reflect which
    studies deposited sequences when, not epidemic change.
    """
    if len(numerators) != len(denominators):
        raise ValueError("numerators and denominators must align")
    if len(numerators) < 3:
        raise ValueError("trend test needs at least three ordered strata")

    scores = list(scores) if scores is not None else [float(i) for i in range(len(numerators))]
    n_total = sum(denominators)
    r_total = sum(numerators)
    if n_total == 0 or r_total in (0, n_total):
        return 0.0, 1.0

    p_bar = r_total / n_total
    s_bar = sum(s * n for s, n in zip(scores, denominators, strict=True)) / n_total
    numer = sum(
        s * (r - n * p_bar) for s, r, n in zip(scores, numerators, denominators, strict=True)
    )
    var = (
        p_bar
        * (1 - p_bar)
        * sum(n * (s - s_bar) ** 2 for s, n in zip(scores, denominators, strict=True))
    )
    if var <= 0:
        return 0.0, 1.0
    z = numer / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2))  # erfc: 1 - erf underflows to 0 for large |z|
    return z, p


def denominator_for_class(
    calls: Iterable[ClassCoverage], drug_class: str, min_coverage: float
) -> int:
    """Count sequences eligible to contribute to ``drug_class``.

    A sequence with no integrase coverage contributes to no INSTI denominator --
    using one global n would turn missing data into apparent susceptibility.
    """
    return sum(1 for call in calls if call.covers_class(drug_class, min_coverage))
