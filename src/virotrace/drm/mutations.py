"""Amino-acid mutation representation and gene-local numbering.

Resistance mutations are reported in GENE-LOCAL amino-acid numbering (RT M184V is
RT codon 184), which is the convention every resistance database uses. HXB2
nucleotide coordinates are used internally for alignment work only. Mixing the two
is a classic HIV-bioinformatics defect, so the boundary is enforced here and
asserted in tests.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Literal

from ..errors import ContractError

Gene = Literal["PR", "RT", "IN"]

#: Protein lengths in amino acids. Asserted against the pinned reference annotation
#: in tests/contract/test_reference.py rather than trusted from this file.
GENE_PROTEIN_LENGTH: dict[str, int] = {"PR": 99, "RT": 560, "IN": 288}

#: Drug classes each gene can inform.
GENE_DRUG_CLASSES: dict[str, tuple[str, ...]] = {
    "PR": ("PI",),
    "RT": ("NRTI", "NNRTI"),
    "IN": ("INSTI",),
}

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
_LABEL_RE = re.compile(r"^(?P<gene>PR|RT|IN):(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z]+|-|\*|\.)$")

CallStatus = Literal["called", "no_call"]


@dataclass(frozen=True, slots=True)
class Mutation:
    """A single amino-acid substitution in gene-local numbering."""

    gene: str
    position: int
    ref_aa: str
    alt_aa: frozenset[str]

    def __post_init__(self) -> None:
        if self.gene not in GENE_PROTEIN_LENGTH:
            raise ContractError(f"unknown gene {self.gene!r}")
        limit = GENE_PROTEIN_LENGTH[self.gene]
        if not 1 <= self.position <= limit:
            raise ContractError(
                f"{self.gene} position {self.position} outside 1..{limit} -- this usually "
                "means an HXB2 nucleotide coordinate leaked into gene-local numbering"
            )
        if not self.alt_aa:
            raise ContractError(f"{self.gene}:{self.position} has no alternate residue")

    @property
    def is_mixture(self) -> bool:
        return len(self.alt_aa) > 1

    @property
    def label(self) -> str:
        alts = "".join(sorted(self.alt_aa))
        return f"{self.gene}:{self.ref_aa}{self.position}{alts}"

    @property
    def position_label(self) -> str:
        return f"{self.gene}:{self.ref_aa}{self.position}"

    @property
    def drug_classes(self) -> tuple[str, ...]:
        return GENE_DRUG_CLASSES[self.gene]


def parse_mutation(label: str) -> Mutation:
    """Parse ``RT:M184V`` or a mixture ``RT:M184IV`` into a :class:`Mutation`.

    Round-trips with :attr:`Mutation.label`; this is property-tested.
    """
    match = _LABEL_RE.match(label.strip().upper())
    if not match:
        raise ContractError(f"unparseable mutation label: {label!r}")
    alt_raw = match["alt"]
    alts = frozenset(alt_raw) if alt_raw not in {"-", "*", "."} else frozenset({alt_raw})
    return Mutation(
        gene=match["gene"],
        position=int(match["pos"]),
        ref_aa=match["ref"],
        alt_aa=alts,
    )


@dataclass(slots=True)
class MutationCalls:
    """The result of calling mutations on one sequence.

    ``covered`` records, per gene, which positions were actually sequenced. An
    uncovered position must resolve to ``no_call`` -- never to wild type. Reporting
    uncovered positions as wild type is the silent false negative that makes a
    resistance survey wrong without ever raising an error.
    """

    strain: str
    mutations: list[Mutation] = field(default_factory=list)
    covered: dict[str, set[int]] = field(default_factory=dict)

    def status_at(self, gene: str, position: int) -> CallStatus:
        if gene not in GENE_PROTEIN_LENGTH:
            raise ContractError(f"unknown gene {gene!r}")
        return "called" if position in self.covered.get(gene, set()) else "no_call"

    def coverage_fraction(self, gene: str) -> float:
        if gene not in GENE_PROTEIN_LENGTH:
            raise ContractError(f"unknown gene {gene!r}")
        return len(self.covered.get(gene, set())) / GENE_PROTEIN_LENGTH[gene]

    def covers_class(self, drug_class: str, min_fraction: float) -> bool:
        """Whether this sequence can contribute to a denominator for ``drug_class``."""
        genes = [g for g, classes in GENE_DRUG_CLASSES.items() if drug_class in classes]
        if not genes:
            raise ContractError(f"unknown drug class {drug_class!r}")
        return any(self.coverage_fraction(g) >= min_fraction for g in genes)

    def labels(self) -> set[str]:
        return {m.label for m in self.mutations}


def collapse_mixture(residues: set[str], policy: str) -> frozenset[str]:
    """Apply the configured mixture policy to an ambiguous codon's residue set.

    ``any``      -- keep every residue; a position is resistant if any call is a DRM.
    ``strict``   -- keep nothing unless the call is unambiguous.
    ``majority`` -- not inferable from sequence alone; treated as ``strict`` and
                    flagged, rather than silently guessing a dominant variant.
    """
    clean = {r for r in residues if r in AMINO_ACIDS}
    if not clean:
        raise ContractError("mixture contains no standard amino acid")
    if policy == "any":
        return frozenset(clean)
    if policy in {"strict", "majority"}:
        return frozenset(clean) if len(clean) == 1 else frozenset()
    raise ContractError(f"unknown mixture policy {policy!r}")
