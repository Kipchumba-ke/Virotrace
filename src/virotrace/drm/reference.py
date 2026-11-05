"""HXB2 reference proteins, derived from the pinned GenBank record.

The gene coordinates live in config; this module only translates what those
coordinates cut out of ``config/reference/hxb2.gb``. It exists so that both the
contract tests and the SDRM importer check against the same translation, instead of
trusting a hand-typed protein sequence.
"""

from __future__ import annotations

from pathlib import Path

from ..acquire.genbank_parse import parse_genbank
from ..errors import ContractError

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


def translate(nucleotides: str) -> str:
    """Translate in frame; a trailing partial codon is a contract error."""
    nucleotides = nucleotides.upper()
    if len(nucleotides) % 3:
        raise ContractError(f"cannot translate {len(nucleotides)} nt: not a multiple of 3")
    return "".join(
        CODON_TABLE.get(nucleotides[i : i + 3], "X") for i in range(0, len(nucleotides), 3)
    )


def load_reference_sequence(path: Path | str) -> str:
    path = Path(path)
    if not path.exists():
        raise ContractError(f"reference record not found: {path}")
    records = list(parse_genbank(path.read_text()))
    if len(records) != 1:
        raise ContractError(f"{path}: expected exactly one record, found {len(records)}")
    return records[0].sequence


def gene_proteins(sequence: str, spans: dict[str, tuple[int, int]]) -> dict[str, str]:
    """Protein per gene from 1-based inclusive nucleotide spans, e.g. ``{"PR": (2253, 2549)}``."""
    return {gene: translate(sequence[start - 1 : end]) for gene, (start, end) in spans.items()}
