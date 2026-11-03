"""The temporal-signal gate.

A time tree from data without temporal signal looks exactly as convincing as one
from data with it, so this decision is automated rather than left to judgement.
"""

from __future__ import annotations

import pytest

from virotrace.phylo.temporal import evaluate_gate, linear_regression


def test_regression_recovers_a_known_slope():
    x = [2010, 2012, 2014, 2016, 2018]
    rate = 1.5e-3
    y = [(xi - 2010) * rate for xi in x]
    slope, _intercept, r2 = linear_regression(x, y)
    assert slope == pytest.approx(rate, rel=1e-9)
    assert r2 == pytest.approx(1.0)


def test_gate_passes_with_clean_temporal_signal():
    x = [2008, 2010, 2012, 2014, 2016, 2018, 2020]
    y = [(xi - 2008) * 1.8e-3 for xi in x]
    signal = evaluate_gate(x, y, min_r2=0.10)
    assert signal.passed is True
    assert signal.implied_rate > 0
    assert "adequate" in signal.reason


def test_gate_fails_on_flat_divergence_and_says_why():
    x = [2008, 2010, 2012, 2014, 2016, 2018]
    y = [0.05, 0.051, 0.049, 0.05, 0.0505, 0.0495]
    signal = evaluate_gate(x, y, min_r2=0.10)
    assert signal.passed is False
    assert "R^2" in signal.reason
    assert "divergence tree" in signal.reason


def test_gate_fails_on_negative_slope():
    x = [2008, 2010, 2012, 2014, 2016]
    y = [0.10, 0.08, 0.06, 0.04, 0.02]
    signal = evaluate_gate(x, y, min_r2=0.0, require_positive_slope=True)
    assert signal.passed is False
    assert "not positive" in signal.reason


def test_identical_dates_are_rejected_rather_than_regressed():
    with pytest.raises(ValueError, match="no temporal axis"):
        evaluate_gate([2015, 2015, 2015], [0.01, 0.02, 0.03])


def test_too_few_points():
    with pytest.raises(ValueError, match="three points"):
        linear_regression([2015, 2016], [0.1, 0.2])
