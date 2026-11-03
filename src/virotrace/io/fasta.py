"""Minimal FASTA I/O with an enforced header policy.

Header policy: a record ID is the ``strain`` and nothing else. Metadata lives in the
TSV, never in the FASTA header -- embedded metadata is how a pipeline ends up with
two disagreeing copies of a sequence's country.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from pathlib import Path

from ..errors import ContractError

STRAIN_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")


def read_fasta(path: Path | str) -> Iterator[tuple[str, str]]:
    """Yield ``(strain, sequence)`` pairs, uppercased, newlines stripped."""
    path = Path(path)
    strain: str | None = None
    chunks: list[str] = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if strain is not None:
                    yield strain, "".join(chunks)
                strain = line[1:].split()[0] if len(line) > 1 else ""
                if not STRAIN_RE.match(strain):
                    raise ContractError(
                        f"{path}: header {strain!r} violates the strain-only header policy"
                    )
                chunks = []
            else:
                if strain is None:
                    raise ContractError(f"{path}: sequence data before any header")
                chunks.append(line.upper())

    if strain is not None:
        yield strain, "".join(chunks)


def write_fasta(path: Path | str, records: Iterable[tuple[str, str]], wrap: int = 70) -> int:
    """Write records, returning the count. Strains are validated and must be unique."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    count = 0

    with path.open("w") as handle:
        for strain, sequence in records:
            if not STRAIN_RE.match(strain):
                raise ContractError(f"invalid strain id for FASTA output: {strain!r}")
            if strain in seen:
                raise ContractError(f"duplicate strain in FASTA output: {strain!r}")
            seen.add(strain)
            handle.write(f">{strain}\n")
            if wrap and wrap > 0:
                for i in range(0, len(sequence), wrap):
                    handle.write(sequence[i : i + wrap] + "\n")
            else:
                handle.write(sequence + "\n")
            count += 1
    return count


def assert_equal_length(records: list[tuple[str, str]]) -> int:
    """Assert an alignment is rectangular; return the alignment length."""
    if not records:
        raise ContractError("alignment is empty")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ContractError(
            f"alignment is ragged: {len(lengths)} distinct lengths {sorted(lengths)}"
        )
    return lengths.pop()
