"""Drug class of an SDRM position.

The SDRM tables list (gene, position, residues) and carry no class, but prevalence is
reported per class and RT holds two of them. The mapping is explicit and closed: an RT
position on a list that is in neither set is an error, not a silent NRTI.
"""

from __future__ import annotations

from ..errors import ContractError

DRUG_CLASSES = ("NRTI", "NNRTI", "PI", "INSTI")

#: RT codons that confer NRTI resistance on the WHO 2009 and Stanford SDRM lists.
NRTI_POSITIONS = frozenset({41, 65, 67, 69, 70, 74, 75, 77, 115, 116, 151, 184, 210, 215, 219})
#: RT codons that confer NNRTI resistance on the same lists.
NNRTI_POSITIONS = frozenset({100, 101, 103, 106, 179, 181, 188, 190, 221, 225, 227, 230, 234, 318})

_GENE_CLASS = {"PR": "PI", "IN": "INSTI"}


def drug_class_of(gene: str, position: int) -> str:
    if gene in _GENE_CLASS:
        return _GENE_CLASS[gene]
    if gene == "RT":
        if position in NRTI_POSITIONS:
            return "NRTI"
        if position in NNRTI_POSITIONS:
            return "NNRTI"
        raise ContractError(f"RT position {position} has no drug class; add it to drm/classes.py")
    raise ContractError(f"unknown gene {gene!r}")


def gene_of_class(drug_class: str) -> str:
    return {"PI": "PR", "NRTI": "RT", "NNRTI": "RT", "INSTI": "IN"}[drug_class]
