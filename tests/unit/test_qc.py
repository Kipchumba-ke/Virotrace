"""QC rules and the record-conservation invariant."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from virotrace.curate.qc import (
    ambiguity_fraction,
    apply_qc,
    gene_coverage,
    has_internal_stop,
    is_frameshifted,
)

GENE_NT = 99  # a short synthetic "gene" for fixtures


def seq(pattern: str, length: int = GENE_NT) -> str:
    return (pattern * (length // len(pattern) + 1))[:length]


def ambiguous_but_covered(n_ambiguous: int = 7) -> str:
    """Enough Ns to trip QC-2 while still passing the QC-1 coverage floor.

    Exercises rule ordering: a sequence that is merely ambiguous must be reported
    as QC-2, not swept up by QC-1.
    """
    bases = list(seq("ACG"))
    for i in range(n_ambiguous):
        bases[i * 5] = "N"
    return "".join(bases)


def test_ambiguity_fraction_ignores_gaps():
    assert ambiguity_fraction("ACGT") == 0.0
    assert ambiguity_fraction("ACGN") == pytest.approx(0.25)
    assert ambiguity_fraction("AC--GN") == pytest.approx(0.25)
    assert ambiguity_fraction("----") == 1.0


def test_gene_coverage_counts_only_unambiguous_bases():
    assert gene_coverage("ACGT" * 25, 100) == pytest.approx(1.0)
    assert gene_coverage("ACGT" * 10, 100) == pytest.approx(0.4)
    assert gene_coverage("NNNN" * 25, 100) == 0.0


def test_gene_coverage_rejects_bad_length():
    with pytest.raises(ValueError):
        gene_coverage("ACGT", 0)


def test_internal_stop_detected_terminal_ignored():
    assert has_internal_stop("ATGTAAATG") is True
    assert has_internal_stop("ATGATGTAA") is False  # terminal stop is fine
    assert has_internal_stop("ATGATGATG") is False


def test_frameshift_detection():
    assert is_frameshifted("ACG---ACG") is False  # 3bp gap, in frame
    assert is_frameshifted("ACG--ACG") is True  # 2bp gap
    assert is_frameshifted("ACGACG") is False


# --------------------------------------------------------------------- apply_qc
def test_records_are_conserved_and_each_exclusion_has_one_reason():
    records = [
        ("good_1", seq("ACG")),
        ("short", "ACGACG"),  # QC-1
        ("ambiguous", ambiguous_but_covered()),  # QC-2
        ("stopped", "ATGTAA" + seq("ACG", 93)),  # QC-3
        ("lab_HXB2", seq("ACG")),  # QC-7
        ("dup_of_good", seq("ACG")),  # QC-6
    ]
    result = apply_qc(
        records,
        gene_length_nt=GENE_NT,
        min_gene_coverage=0.80,
        max_ambiguity_frac=0.05,
        lab_strains=frozenset({"lab_HXB2"}),
    )

    result.assert_conserved()
    assert result.n_in == len(records)
    assert result.n_in == result.n_out + result.n_excluded

    reasons = {entry.strain: entry.rule for entry in result.ledger}
    assert reasons["short"] == "QC-1"
    assert reasons["ambiguous"] == "QC-2"
    assert reasons["stopped"] == "QC-3"
    assert reasons["lab_HXB2"] == "QC-7"
    assert reasons["dup_of_good"] == "QC-6"
    assert result.kept == ["good_1"]

    # One primary reason per record, so the QC waterfall cannot double-count.
    assert len(reasons) == len(result.ledger)


def test_ledger_records_observed_value_and_threshold():
    result = apply_qc(
        [("short", "ACGACG")],
        gene_length_nt=GENE_NT,
        min_gene_coverage=0.80,
        max_ambiguity_frac=0.05,
    )
    entry = result.ledger[0]
    assert entry.rule == "QC-1"
    assert entry.threshold == 0.80
    assert isinstance(entry.observed, float) and entry.observed < 0.80


@given(n_good=st.integers(0, 20), n_short=st.integers(0, 20))
def test_conservation_is_a_property(n_good, n_short):
    records = [(f"g{i}", seq("ACG")) for i in range(n_good)]
    records += [(f"s{i}", "AC") for i in range(n_short)]
    # Distinguish the "good" duplicates so dedup does not confound the count.
    records = [
        (name, s if name.startswith("s") else s[: GENE_NT - i] + "ACG" * 1)
        for i, (name, s) in enumerate(records)
    ]
    result = apply_qc(
        records,
        gene_length_nt=GENE_NT,
        min_gene_coverage=0.0,
        max_ambiguity_frac=1.0,
        exclude_internal_stop=False,
        exclude_frameshift=False,
    )
    assert result.n_in == result.n_out + result.n_excluded
