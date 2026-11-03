"""Temporal signal gate.

A time-calibrated tree built from data with no temporal signal looks exactly as
convincing as one built from data with it. That is why this gate is automated rather
than left to judgement at writing time: if root-to-tip divergence does not increase
with sampling date, the pipeline publishes the divergence tree and states why.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TemporalSignal:
    slope: float
    intercept: float
    r_squared: float
    n: int
    passed: bool
    reason: str

    @property
    def implied_rate(self) -> float:
        """Slope in substitutions/site/year, i.e. a crude clock rate estimate."""
        return self.slope


def linear_regression(x: Sequence[float], y: Sequence[float]) -> tuple[float, float, float]:
    """Ordinary least squares. Returns ``(slope, intercept, r_squared)``."""
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must be the same length")
    if n < 3:
        raise ValueError("need at least three points for a regression")

    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sxx = sum((xi - mean_x) ** 2 for xi in x)
    sxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y, strict=True))
    syy = sum((yi - mean_y) ** 2 for yi in y)

    if sxx == 0:
        raise ValueError("all sampling dates are identical -- no temporal axis")

    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    r_squared = 0.0 if syy == 0 else (sxy**2) / (sxx * syy)
    return slope, intercept, r_squared


def evaluate_gate(
    dates: Sequence[float],
    root_to_tip: Sequence[float],
    *,
    min_r2: float = 0.10,
    require_positive_slope: bool = True,
) -> TemporalSignal:
    """Decide whether a time-calibrated tree is defensible for this dataset."""
    slope, intercept, r2 = linear_regression(dates, root_to_tip)

    reasons: list[str] = []
    if require_positive_slope and slope <= 0:
        reasons.append(
            f"root-to-tip slope is {slope:.3e} (not positive): divergence does not "
            "increase with sampling date"
        )
    if r2 < min_r2:
        reasons.append(f"root-to-tip R^2 is {r2:.3f}, below the required {min_r2:.3f}")

    passed = not reasons
    return TemporalSignal(
        slope=slope,
        intercept=intercept,
        r_squared=r2,
        n=len(dates),
        passed=passed,
        reason=(
            "temporal signal adequate"
            if passed
            else "; ".join(reasons) + " -- emitting divergence tree instead of a time tree"
        ),
    )
