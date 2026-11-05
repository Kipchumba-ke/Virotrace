"""Project one aligned sequence onto the reference's own coordinates.

A pairwise-with-reference alignment has extra columns wherever the sequence has bases
the reference lacks. For resistance calling we want every sequence in the same
rectangular coordinate system (one column per HXB2 nucleotide), so those insertion
columns are removed -- but not forgotten. An insertion whose length is not a multiple
of three is a frameshift, and a plain ``--keeplength`` alignment would hide it, so the
insertions are returned alongside the projected sequence and QC reads them.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..errors import ContractError

GAP_CHARS = frozenset("-.")


@dataclass(frozen=True, slots=True)
class Insertion:
    """``length`` extra bases in the sequence, after reference position ``after_ref_pos``."""

    after_ref_pos: int
    length: int


@dataclass(frozen=True, slots=True)
class Projection:
    #: one character per reference position in [start, end]; '-' where the sequence has
    #: no base (a deletion, or the sequence simply does not extend that far)
    projected: str
    insertions: tuple[Insertion, ...]


def project_to_reference(
    aligned_reference: str, aligned_sequence: str, *, start: int, end: int
) -> Projection:
    """Project onto reference positions ``start..end`` (1-based, inclusive).

    ``aligned_reference`` and ``aligned_sequence`` are rows of the same pairwise
    alignment. Reference positions are counted from the first base of the reference
    row, so ``start``/``end`` are positions within that row's ungapped sequence.
    """
    if len(aligned_reference) != len(aligned_sequence):
        raise ContractError(
            f"aligned rows differ in length ({len(aligned_reference)} vs {len(aligned_sequence)})"
        )
    if not 1 <= start <= end:
        raise ContractError(f"invalid reference span {start}..{end}")

    width = end - start + 1
    projected = ["-"] * width
    insertions: list[Insertion] = []
    ref_pos = 0
    run = 0  # length of the insertion run currently being read

    for ref_char, seq_char in zip(aligned_reference, aligned_sequence, strict=True):
        seq_base = "-" if seq_char in GAP_CHARS else seq_char.upper()
        if ref_char in GAP_CHARS:
            if seq_base != "-":
                run += 1
            continue
        if run:
            # Only insertions strictly inside the span matter: bases hanging off either
            # end are flanking sequence, not a defect in the target region.
            if start <= ref_pos < end:
                insertions.append(Insertion(ref_pos, run))
            run = 0
        ref_pos += 1
        if start <= ref_pos <= end:
            projected[ref_pos - start] = seq_base

    if ref_pos < end:
        raise ContractError(
            f"reference row has only {ref_pos} bases, but the span runs to position {end}"
        )
    return Projection("".join(projected), tuple(insertions))
