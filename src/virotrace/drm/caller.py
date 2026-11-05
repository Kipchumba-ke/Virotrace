"""Codon-aware mutation calling from a reference-projected alignment.

Input is a sequence already projected onto HXB2 pol coordinates (see
``align.project``). For each gene this reads the sequence codon by codon and decides,
for every codon position, one of: called (one residue or a mixture), deleted, or
no-call. The rule that matters most is the last one: a position the sequence does not
cover is a no-call and is NEVER reported as wild type, because that silent false
negative would understate resistance without raising any error.

Two defects disqualify a whole gene rather than a position, since every amino-acid
call from it would be unreliable: an indel whose length is not a multiple of three
(frameshift) and an in-frame internal stop codon. A disqualified gene contributes no
covered positions; the QC stage records why.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from ..align.project import Insertion
from ..errors import ContractError
from .mutations import GENE_PROTEIN_LENGTH, Mutation, MutationCalls, collapse_mixture
from .reference import CODON_TABLE

IUPAC: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}  # fmt: skip

CodonKind = str  # "aa" | "stop" | "too_ambiguous" | "invalid"


def classify_codon(codon: str, max_residues: int) -> tuple[CodonKind, frozenset[str]]:
    """Residues a (possibly ambiguous) codon can encode.

    A codon that could encode more than ``max_residues`` amino acids carries no usable
    information (``NNN`` is the extreme) and is reported as ``too_ambiguous``. A codon
    whose every reading is a stop is a stop; a mixture that merely *includes* a stop
    reading is judged on its amino-acid readings.
    """
    codon = codon.upper()
    if len(codon) != 3 or any(base not in IUPAC for base in codon):
        return "invalid", frozenset()
    first, second, third = (IUPAC[base] for base in codon)
    residues = {CODON_TABLE[a + b + c] for a in first for b in second for c in third}
    residues.discard("*")
    if not residues:
        return "stop", frozenset()
    if len(residues) > max_residues:
        return "too_ambiguous", frozenset()
    return "aa", frozenset(residues)


@dataclass(slots=True)
class GeneResult:
    gene: str
    #: gene-local position -> residues called there (more than one means a mixture)
    called: dict[int, frozenset[str]] = field(default_factory=dict)
    deleted: set[int] = field(default_factory=set)
    internal_stops: list[int] = field(default_factory=list)
    frameshift: bool = False
    n_too_ambiguous: int = 0

    @property
    def usable(self) -> bool:
        return not self.frameshift and not self.internal_stops

    @property
    def covered(self) -> set[int]:
        return (set(self.called) | self.deleted) if self.usable else set()


@dataclass(slots=True)
class SequenceCalls:
    strain: str
    genes: dict[str, GeneResult]
    calls: MutationCalls


def _gap_runs_inside(projected: str, first: int, last: int) -> list[tuple[int, int]]:
    """``(start_index, length)`` of '-' runs strictly between the first and last base.

    Measured over the whole projected sequence, so a deletion that straddles a gene
    boundary keeps its true length instead of being cut into two misleading pieces.
    """
    runs: list[tuple[int, int]] = []
    current = 0
    for index in range(first, last + 1):
        if projected[index] == "-":
            current += 1
        elif current:
            runs.append((index - current, current))
            current = 0
    return runs


def insertions_by_gene(
    insertions: Iterable[Insertion],
    gene_spans: Mapping[str, tuple[int, int]],
    pol_start: int,
) -> dict[str, list[int]]:
    """Group insertion lengths by the gene they fall in.

    ``Insertion.after_ref_pos`` counts from the first base of the pol reference, so it is
    converted to an HXB2 coordinate before comparing with the gene spans.
    """
    out: dict[str, list[int]] = {gene: [] for gene in gene_spans}
    for insertion in insertions:
        hxb2 = pol_start - 1 + insertion.after_ref_pos
        for gene, (g_start, g_end) in gene_spans.items():
            if g_start <= hxb2 < g_end:
                out[gene].append(insertion.length)
    return out


def call_projected(
    strain: str,
    projected: str,
    *,
    pol_start: int,
    gene_spans: Mapping[str, tuple[int, int]],
    reference_protein: Mapping[str, str],
    insertion_lengths: Mapping[str, list[int]] | None = None,
    mixture_policy: str = "any",
    max_ambiguous_residues: int = 4,
) -> SequenceCalls:
    """Call mutations for one projected sequence.

    ``insertion_lengths`` maps gene -> lengths of insertions the sequence carries inside
    that gene (from ``Projection.insertions``); any length not divisible by three is a
    frameshift.
    """
    bases = [i for i, char in enumerate(projected) if char != "-"]
    first, last = (bases[0], bases[-1]) if bases else (0, -1)
    runs = _gap_runs_inside(projected, first, last) if bases else []
    genes: dict[str, GeneResult] = {}

    for gene, (g_start, g_end) in gene_spans.items():
        if gene not in GENE_PROTEIN_LENGTH:
            raise ContractError(f"unknown gene {gene!r}")
        n_codons = GENE_PROTEIN_LENGTH[gene]
        if g_end - g_start + 1 != 3 * n_codons:
            raise ContractError(f"{gene} span {g_start}..{g_end} does not match {n_codons} codons")
        reference = reference_protein[gene]
        result = GeneResult(gene)

        offset = g_start - pol_start
        gene_slice = projected[offset : offset + 3 * n_codons]
        if len(gene_slice) != 3 * n_codons:
            raise ContractError(f"{strain}: projected sequence too short for {gene}")

        # Frame: gaps strictly inside the sequence's extent, and insertions, must be
        # whole codons. Leading/trailing gaps are just missing coverage.
        gene_lo, gene_hi = offset, offset + 3 * n_codons - 1
        result.frameshift = any(
            length % 3 and start <= gene_hi and start + length - 1 >= gene_lo
            for start, length in runs
        )
        if any(length % 3 for length in (insertion_lengths or {}).get(gene, [])):
            result.frameshift = True

        for index in range(n_codons):
            codon = gene_slice[3 * index : 3 * index + 3]
            position = index + 1
            col0 = offset + 3 * index
            if col0 + 2 < first or col0 > last:
                continue  # the sequence does not reach this codon: no call
            if codon == "---":
                result.deleted.add(position)
                continue
            if "-" in codon:
                continue  # partly missing (edge of coverage, or part of a frameshift)
            kind, residues = classify_codon(codon, max_ambiguous_residues)
            if kind == "stop":
                if reference[index] != "*":
                    result.internal_stops.append(position)
            elif kind == "too_ambiguous":
                result.n_too_ambiguous += 1
            elif kind == "aa":
                result.called[position] = residues
        genes[gene] = result

    return SequenceCalls(
        strain, genes, _to_mutation_calls(strain, genes, reference_protein, mixture_policy)
    )


def _to_mutation_calls(
    strain: str,
    genes: Mapping[str, GeneResult],
    reference_protein: Mapping[str, str],
    mixture_policy: str,
) -> MutationCalls:
    calls = MutationCalls(strain=strain)
    for gene, result in genes.items():
        covered = result.covered
        if not covered:
            calls.covered[gene] = set()
            continue
        covered = set(covered)
        for position, residues in sorted(result.called.items()):
            reference = reference_protein[gene][position - 1]
            kept = collapse_mixture(set(residues), mixture_policy)
            if not kept:
                # strict policy and an ambiguous codon: not a usable call, so not covered
                covered.discard(position)
                continue
            if kept != {reference}:
                calls.mutations.append(Mutation(gene, position, reference, kept))
        calls.covered[gene] = covered
    return calls
