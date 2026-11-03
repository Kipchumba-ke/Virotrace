"""Wilson intervals, per-class denominators, suppression, and trend."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from virotrace.drm.prevalence import (
    SUPPRESSED,
    cochran_armitage_trend,
    estimate_prevalence,
    wilson_ci,
)


def test_wilson_interval_contains_the_estimate():
    low, high = wilson_ci(12, 100)
    assert low < 0.12 < high


def test_wilson_handles_zero_successes_without_a_degenerate_interval():
    low, high = wilson_ci(0, 30)
    assert low == 0.0
    assert 0 < high < 0.2  # Wald would give [0, 0], which is the bug this avoids


def test_wilson_handles_all_successes():
    """At p=1 the upper bound is analytically exactly 1, not 1 minus a few ulps."""
    low, high = wilson_ci(30, 30)
    assert high == 1.0
    assert 0.8 < low < 1.0


@given(
    total=st.integers(min_value=1, max_value=5000),
    frac=st.floats(min_value=0.0, max_value=1.0),
)
def test_wilson_bounds_are_valid_and_ordered(total, frac):
    successes = round(frac * total)
    low, high = wilson_ci(successes, total)
    assert 0.0 <= low <= high <= 1.0
    assert low <= successes / total <= high


def test_wilson_rejects_impossible_counts():
    with pytest.raises(ValueError):
        wilson_ci(5, 3)
    with pytest.raises(ValueError):
        wilson_ci(1, 0)


# ------------------------------------------------------------------- suppression
def test_small_denominator_is_suppressed_not_plotted():
    """An unsuppressed 100% from n=2 is the classic dashboard embarrassment."""
    row = estimate_prevalence(numerator=2, denominator=2, drug_class="INSTI", min_stratum_n=10)
    assert row.status == SUPPRESSED
    assert row.proportion is None
    assert row.ci_low is None
    # The stratum is still reported, so the reader sees it existed.
    assert row.denominator == 2


def test_adequate_denominator_is_estimated():
    row = estimate_prevalence(numerator=5, denominator=50, drug_class="NNRTI", min_stratum_n=10)
    assert row.status == "estimated"
    assert row.proportion == pytest.approx(0.1)
    assert row.ci_low is not None and row.ci_high is not None
    assert row.ci_low < row.proportion < row.ci_high


def test_numerator_above_denominator_is_an_error():
    """Guards the per-class denominator logic, where this is the likely failure."""
    with pytest.raises(ValueError, match="exceeds denominator"):
        estimate_prevalence(numerator=11, denominator=10, drug_class="NRTI")


# ------------------------------------------------------------------------ trend
def test_trend_detects_a_monotone_increase():
    z, p = cochran_armitage_trend([2, 8, 20], [100, 100, 100])
    assert z > 0
    assert p < 0.01


def test_trend_is_flat_for_constant_proportions():
    z, p = cochran_armitage_trend([10, 10, 10], [100, 100, 100])
    assert abs(z) < 1e-9
    assert p == pytest.approx(1.0, abs=1e-9)


def test_trend_needs_three_strata():
    with pytest.raises(ValueError, match="three ordered strata"):
        cochran_armitage_trend([1, 2], [10, 10])
