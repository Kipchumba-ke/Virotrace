"""WHO surveillance drug resistance mutation (SDRM) tagging.

Two mutation sets serve different purposes and must not be substituted for each
other:

* the full HIVdb list with penalty scores -- per-sequence clinical-style
  interpretation, sensitive to subtype polymorphism;
* the WHO SDRM list -- curated to be comparable across populations and to exclude
  common polymorphisms, which makes it the correct set for a *surveillance*
  prevalence claim.

Reporting HIVdb penalty scores as a prevalence metric would not be comparable to
published regional estimates. SDRM prevalence is the headline; HIVdb interpretation
is the secondary per-sequence view.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ..errors import ContractError
from .mutations import GENE_PROTEIN_LENGTH, Mutation


@dataclass(frozen=True, slots=True)
class SDRMList:
    version: str
    #: (gene, position) -> set of resistance-conferring residues
    entries: dict[tuple[str, int], frozenset[str]]

    def is_sdrm(self, mutation: Mutation) -> bool:
        residues = self.entries.get((mutation.gene, mutation.position))
        if residues is None:
            return False
        # Mixture policy 'any': resistant if any called residue is on the list.
        return bool(mutation.alt_aa & residues)

    def matched_residues(self, mutation: Mutation) -> frozenset[str]:
        residues = self.entries.get((mutation.gene, mutation.position), frozenset())
        return frozenset(mutation.alt_aa & residues)

    def __len__(self) -> int:
        return len(self.entries)


def load_sdrm_list(path: Path | str, version: str) -> SDRMList:
    """Load a TSV with columns ``gene``, ``position``, ``residues``.

    ``residues`` is a string of one-letter codes, e.g. ``IV`` for positions where
    either isoleucine or valine is an SDRM.
    """
    path = Path(path)
    if not path.exists():
        raise ContractError(f"SDRM list not found: {path}")

    entries: dict[tuple[str, int], frozenset[str]] = {}
    with path.open(newline="") as handle:
        # csv.DictReader has no comment support, and a '#' header line would
        # otherwise be read as a data row with gene='# ...'.
        rows = (line for line in handle if not line.lstrip().startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        required = {"gene", "position", "residues"}
        if not required <= set(reader.fieldnames or []):
            raise ContractError(
                f"SDRM list {path} must have columns {sorted(required)}, found {reader.fieldnames}"
            )
        for row in reader:
            gene = row["gene"].strip().upper()
            if gene not in GENE_PROTEIN_LENGTH:
                raise ContractError(f"SDRM list {path}: unknown gene {gene!r}")
            position = int(row["position"])
            if not 1 <= position <= GENE_PROTEIN_LENGTH[gene]:
                raise ContractError(
                    f"SDRM list {path}: {gene} position {position} outside "
                    f"1..{GENE_PROTEIN_LENGTH[gene]}"
                )
            residues = frozenset(row["residues"].strip().upper())
            entries[(gene, position)] = entries.get((gene, position), frozenset()) | residues

    if not entries:
        raise ContractError(f"SDRM list {path} is empty")
    return SDRMList(version=version, entries=entries)
