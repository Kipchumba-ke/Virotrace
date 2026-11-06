"""Temporal signal gate.

A time-calibrated tree built from data with no temporal signal looks exactly as
convincing as one built from data with it. That is why this gate is automated rather
than left to judgement at writing time: if root-to-tip divergence does not increase
with sampling date, the pipeline publishes the divergence tree and states why.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any


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


# --- root-to-tip with the best root ---------------------------------------------------------


def decimal_year(date: str) -> float | None:
    """``2010-05-14`` -> 2010.36; ``2010-05`` -> mid-month; ``2010`` -> mid-year; else None."""
    import calendar
    from datetime import date as _date

    # augur-style ambiguity: 2011-11-XX is a month, 2011-XX-XX a year
    parts = [p for p in date.strip().split("-") if p]
    while parts and parts[-1].upper() in {"XX", "XXXX", "?"}:
        parts.pop()
    try:
        year = int(parts[0])
        if len(parts) >= 3:
            d = _date(year, int(parts[1]), int(parts[2]))
            span = 366 if calendar.isleap(year) else 365
            return year + (d.timetuple().tm_yday - 0.5) / span
        if len(parts) == 2:
            return year + (int(parts[1]) - 0.5) / 12
        return year + 0.5
    except (ValueError, IndexError):
        return None


def best_root_to_tip(
    tree: Any,
    dates: Sequence[float],
    *,
    subset: Sequence[bool] | None = None,
    grid: int = 9,
    require_positive_slope: bool = True,
) -> tuple[list[float], list[float]]:
    """Root-to-tip distances for the root (any edge, ``grid`` positions each) with the best R^2.

    ``dates`` are in tip order (``tree.tips``). With ``subset``, only those tips enter the
    regression but the root may still sit anywhere in the full tree. Returns
    ``(dates_used, root_to_tip_used)``.
    """
    import numpy as np

    from .newick import tip_distances

    y_all = np.asarray(dates, dtype=float)
    keep = np.ones(len(y_all), dtype=bool) if subset is None else np.asarray(subset, dtype=bool)
    if keep.sum() < 3:
        raise ValueError("need at least three tips for root-to-tip")
    dist, lo, hi = tip_distances(tree)
    y = y_all[keep]
    yc = y - y.mean()
    var_y = float((yc**2).sum())
    if var_y == 0:
        raise ValueError("all sampling dates are identical -- no temporal axis")
    fractions = np.linspace(0.0, 1.0, grid)
    best: tuple[float, np.ndarray] | None = None
    best_any: tuple[float, np.ndarray] | None = None
    n_tips = len(y_all)
    for node in range(1, tree.n_nodes):
        length = tree.length[node]
        p = tree.parent[node]
        inside = np.zeros(n_tips, dtype=bool)
        inside[lo[node] : hi[node]] = True
        base = np.where(inside, dist[node], dist[p] + length)[keep]
        step = np.where(inside, 1.0, -1.0)[keep] * length
        d = base[None, :] + fractions[:, None] * step[None, :]  # one row per root position
        dc = d - d.mean(axis=1, keepdims=True)
        var_d = (dc**2).sum(axis=1)
        cov = dc @ yc
        with np.errstate(divide="ignore", invalid="ignore"):
            r2 = np.where(var_d > 0, cov * cov / (var_d * var_y), -1.0)
        k = int(np.argmax(r2))
        if r2[k] >= 0 and (best_any is None or r2[k] > best_any[0]):
            best_any = (float(r2[k]), d[k])
        pos = np.where(cov > 0, r2, -1.0)
        k = int(np.argmax(pos))
        if pos[k] >= 0 and (best is None or pos[k] > best[0]):
            best = (float(pos[k]), d[k])
    chosen = best if (best is not None or not require_positive_slope) else best_any
    if chosen is None:
        chosen = best_any
    if chosen is None:
        raise ValueError("no valid root found")
    return y.tolist(), chosen[1].tolist()


def permutation_p(
    tree: Any,
    dates: Sequence[float],
    observed_r2: float,
    *,
    subset: Sequence[bool] | None = None,
    n: int = 100,
    seed: int = 0,
    require_positive_slope: bool = True,
) -> float:
    """Date-randomisation test: how often do shuffled dates reach the observed R^2?

    Each shuffle repeats the whole best-root search, so the bias from choosing the root
    that maximises R^2 is in the null as well as in the observed value.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    y = np.asarray(dates, dtype=float)
    keep = np.ones(len(y), dtype=bool) if subset is None else np.asarray(subset, dtype=bool)
    idx = np.flatnonzero(keep)
    hits = 0
    for _ in range(n):
        shuffled = y.copy()
        shuffled[idx] = rng.permutation(y[idx])
        used, d = best_root_to_tip(
            tree, shuffled, subset=keep, require_positive_slope=require_positive_slope
        )
        if linear_regression(used, d)[2] >= observed_r2:
            hits += 1
    return (hits + 1) / (n + 1)
