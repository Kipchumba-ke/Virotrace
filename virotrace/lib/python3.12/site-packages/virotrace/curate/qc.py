"""Sequence quality control (rules QC-1 .. QC-7).

Every excluded record lands in the ledger with the rule that excluded it and the
observed value. No sequence disappears silently -- see docs/04 section 2.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

AMBIGUOUS = set("NRYKMSWBDHVXn-")
UNAMBIGUOUS = set("ACGTacgt")
STOP_CODONS = {"TAA", "TAG", "TGA"}

QC_RULES: dict[str, str] = {
    "QC-1": "insufficient gene coverage",
    "QC-2": "ambiguity fraction above threshold",
    "QC-3": "internal stop codon",
    "QC-4": "frameshift relative to reference",
    "QC-5": "APOBEC hypermutation",
    "QC-6": "duplicate sequence",
    "QC-7": "laboratory strain",
}


@dataclass(frozen=True, slots=True)
class LedgerEntry:
    """One line of the exclusion ledger."""

    strain: str
    rule: str
    reason: str
    observed: float | str | None = None
    threshold: float | str | None = None


@dataclass(slots=True)
class QCResult:
    kept: list[str]
    ledger: list[LedgerEntry]
    n_in: int

    @property
    def n_out(self) -> int:
        return len(self.kept)

    @property
    def n_excluded(self) -> int:
        return len(self.ledger)

    def assert_conserved(self) -> None:
        """n_in == n_out + n_excluded.

        The commonest silent defect in a filtering pipeline is losing records with
        no ledger entry, so this is checked at the stage boundary, not just in tests.
        """
        if self.n_in != self.n_out + self.n_excluded:
            raise AssertionError(
                f"record conservation violated: in={self.n_in} out={self.n_out} "
                f"excluded={self.n_excluded}"
            )


def ambiguity_fraction(sequence: str) -> float:
    """Fraction of non-ACGT characters among ungapped positions."""
    ungapped = [base for base in sequence if base not in {"-", ".", " "}]
    if not ungapped:
        return 1.0
    ambiguous = sum(1 for base in ungapped if base.upper() not in {"A", "C", "G", "T"})
    return ambiguous / len(ungapped)


def gene_coverage(sequence: str, gene_length_nt: int) -> float:
    """Fraction of a gene's nucleotides that are unambiguously sequenced."""
    if gene_length_nt <= 0:
        raise ValueError("gene_length_nt must be positive")
    informative = sum(1 for base in sequence if base.upper() in {"A", "C", "G", "T"})
    return min(informative / gene_length_nt, 1.0)


def has_internal_stop(sequence: str) -> bool:
    """Whether an in-frame internal stop codon is present (terminal stop ignored)."""
    clean = sequence.replace("-", "").upper()
    codons = [clean[i : i + 3] for i in range(0, len(clean) - len(clean) % 3, 3)]
    return any(codon in STOP_CODONS for codon in codons[:-1])


def is_frameshifted(sequence: str) -> bool:
    """Whether an aligned sequence carries an indel that is not a multiple of three."""
    runs: list[int] = []
    current = 0
    for base in sequence:
        if base == "-":
            current += 1
        elif current:
            runs.append(current)
            current = 0
    if current:
        runs.append(current)
    return any(run % 3 for run in runs)


def apply_qc(
    records: Iterable[tuple[str, str]],
    *,
    gene_length_nt: int,
    min_gene_coverage: float,
    max_ambiguity_frac: float,
    exclude_internal_stop: bool = True,
    exclude_frameshift: bool = True,
    lab_strains: frozenset[str] = frozenset(),
) -> QCResult:
    """Run QC rules in order over ``(strain, sequence)`` pairs.

    The first rule a record fails is its primary ``qc_reason``; later rules are not
    evaluated for that record, so the QC waterfall figure cannot double-count.
    """
    records = list(records)
    kept: list[str] = []
    ledger: list[LedgerEntry] = []
    seen: dict[str, str] = {}

    for strain, sequence in records:
        if lab_strains and strain in lab_strains:
            ledger.append(LedgerEntry(strain, "QC-7", QC_RULES["QC-7"], strain))
            continue

        coverage = gene_coverage(sequence, gene_length_nt)
        if coverage < min_gene_coverage:
            ledger.append(
                LedgerEntry(strain, "QC-1", QC_RULES["QC-1"], round(coverage, 4), min_gene_coverage)
            )
            continue

        ambiguity = ambiguity_fraction(sequence)
        if ambiguity > max_ambiguity_frac:
            ledger.append(
                LedgerEntry(
                    strain, "QC-2", QC_RULES["QC-2"], round(ambiguity, 4), max_ambiguity_frac
                )
            )
            continue

        if exclude_internal_stop and has_internal_stop(sequence):
            ledger.append(LedgerEntry(strain, "QC-3", QC_RULES["QC-3"]))
            continue

        if exclude_frameshift and is_frameshifted(sequence):
            ledger.append(LedgerEntry(strain, "QC-4", QC_RULES["QC-4"]))
            continue

        fingerprint = sequence.replace("-", "").upper()
        if fingerprint in seen:
            ledger.append(LedgerEntry(strain, "QC-6", QC_RULES["QC-6"], observed=seen[fingerprint]))
            continue
        seen[fingerprint] = strain
        kept.append(strain)

    result = QCResult(kept=kept, ledger=ledger, n_in=len(records))
    result.assert_conserved()
    return result
