"""MAFFT wrapper: align a chunk of sequences against the reference.

Each chunk is aligned to the reference only (``--addfragments`` onto a one-sequence
alignment), never to the other chunks: resistance calling needs every sequence placed
against HXB2, not against each other. Insertions are kept (no ``--keeplength``) so
frameshifts stay visible; ``project.project_to_reference`` then puts every sequence
into the same rectangular coordinate system.

The process call is injected (as the NCBI transport is), so everything around it is
testable without MAFFT installed.
"""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Callable, Iterable, Iterator, Sequence
from pathlib import Path

from ..errors import ContractError
from ..io.fasta import write_fasta

#: (reference_fasta, fragments_fasta, threads) -> alignment as FASTA text
Runner = Callable[[Path, Path, int], str]


def mafft_runner(reference: Path, fragments: Path, threads: int) -> str:
    command = [
        "mafft",
        "--quiet",
        "--thread",
        str(threads),
        "--6merpair",
        "--addfragments",
        str(fragments),
        str(reference),
    ]
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise ContractError("mafft is not installed or not on PATH") from exc
    if process.returncode != 0:
        raise ContractError(f"mafft exited with {process.returncode}: {process.stderr[-500:]}")
    return process.stdout


def parse_aligned_fasta(text: str) -> dict[str, str]:
    """Aligned FASTA to ``{name: row}``; lines are joined and names are kept verbatim."""
    rows: dict[str, list[str]] = {}
    name: str | None = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in rows:
                raise ContractError(f"duplicate sequence name in alignment: {name}")
            rows[name] = []
        elif name is None:
            raise ContractError("alignment text does not start with a FASTA header")
        else:
            rows[name].append(line)
    return {key: "".join(parts) for key, parts in rows.items()}


def chunked(items: Iterable[tuple[str, str]], size: int) -> Iterator[list[tuple[str, str]]]:
    if size < 1:
        raise ContractError("chunk size must be at least 1")
    chunk: list[tuple[str, str]] = []
    for item in items:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def align_chunk(
    reference: tuple[str, str],
    chunk: Sequence[tuple[str, str]],
    *,
    threads: int = 1,
    runner: Runner = mafft_runner,
) -> tuple[str, dict[str, str]]:
    """Align one chunk. Returns ``(aligned_reference_row, {strain: aligned_row})``.

    Every input strain must come back exactly once, or the run is aborted: a sequence
    that silently vanishes here would also vanish from the prevalence denominator.
    """
    ref_name, ref_seq = reference
    with tempfile.TemporaryDirectory(prefix="virotrace_mafft_") as tmp:
        ref_path, frag_path = Path(tmp) / "reference.fasta", Path(tmp) / "fragments.fasta"
        write_fasta(ref_path, [(ref_name, ref_seq)])
        write_fasta(frag_path, list(chunk))
        text = runner(ref_path, frag_path, threads)

    rows = parse_aligned_fasta(text)
    if ref_name not in rows:
        raise ContractError(f"reference {ref_name!r} missing from the alignment output")
    aligned_ref = rows.pop(ref_name)
    expected = {name for name, _ in chunk}
    if set(rows) != expected:
        missing = sorted(expected - set(rows))[:5]
        extra = sorted(set(rows) - expected)[:5]
        raise ContractError(
            f"alignment output does not match its input: missing {missing}, unexpected {extra}"
        )
    widths = {len(aligned_ref), *(len(row) for row in rows.values())}
    if len(widths) != 1:
        raise ContractError("alignment rows are not all the same length")
    return aligned_ref, rows
