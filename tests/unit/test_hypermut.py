"""APOBEC hypermutation screen: sensitivity AND specificity.

The specificity half is the one people forget. A screen that flags everything
removes real sequences and biases the corpus just as surely as no screen at all.
"""

from __future__ import annotations

import pytest

from virotrace.curate.hypermut import fisher_exact_greater, screen_hypermutation


def build_pair() -> tuple[str, str]:
    """Reference with APOBEC and control contexts; query edited only in context."""
    # GGA repeated: position 0 G is followed by G (APOBEC context).
    # GCT repeated: position 0 G is followed by C (control context).
    reference = ("GGA" * 40) + ("GCT" * 40)
    edited = ("AGA" * 40) + ("GCT" * 40)
    return reference, edited


def test_flags_edited_sequence():
    reference, edited = build_pair()
    result = screen_hypermutation("edited", edited, reference, p_threshold=0.05)
    assert result.flagged is True
    assert result.context_rate > result.control_rate
    assert result.p_value < 0.05


def test_spares_unedited_parent():
    reference, _ = build_pair()
    result = screen_hypermutation("parent", reference, reference, p_threshold=0.05)
    assert result.flagged is False
    assert result.context_mutated == 0


def test_spares_sequence_with_g_to_a_outside_apobec_context():
    """G-to-A alone is not hypermutation; the context ratio is what matters."""
    reference = ("GGA" * 40) + ("GCT" * 40)
    # Edit only the control-context G positions.
    query = ("GGA" * 40) + ("ACT" * 40)
    result = screen_hypermutation("control_edited", query, reference, p_threshold=0.05)
    assert result.flagged is False
    assert result.control_mutated > 0


def test_requires_aligned_lengths():
    with pytest.raises(ValueError, match="aligned"):
        screen_hypermutation("x", "ACGT", "ACG", p_threshold=0.05)


def test_fisher_exact_edges():
    assert fisher_exact_greater(0, 0, 0, 0) == 1.0
    assert fisher_exact_greater(10, 0, 0, 10) < 0.001
    assert fisher_exact_greater(5, 5, 5, 5) == pytest.approx(0.6718, abs=1e-3)
    with pytest.raises(ValueError):
        fisher_exact_greater(-1, 0, 0, 0)
