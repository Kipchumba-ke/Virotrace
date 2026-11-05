"""The align stage with an injected aligner: chunking, projection, conservation."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from virotrace.align.mafft import align_chunk, chunked, parse_aligned_fasta
from virotrace.align.stage import run_align
from virotrace.drm.reference import CODON_TABLE
from virotrace.errors import ContractError
from virotrace.io.fasta import read_fasta, write_fasta
from virotrace.io.tables import read_tsv
from virotrace.settings import load_settings

POL_START, POL_END = 2085, 5096
SENSE = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")


def _hxb2_like_genbank(path: Path) -> str:
    """A one-record GenBank file whose sequence has a known pol region."""
    rng = random.Random(7)
    pol = "".join(rng.choice(SENSE) for _ in range(1003)) + "TAA"
    seq = "A" * (POL_START - 1) + pol + "C" * (9719 - POL_END)
    lines = ["LOCUS       K03455                9719 bp    DNA     linear   VRL 01-JAN-2000"]
    lines += ["VERSION     K03455.1", "FEATURES             Location/Qualifiers"]
    lines += ["     source          1..9719", "ORIGIN"]
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60].lower()
        lines.append(f"{i + 1:>9} " + " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10)))
    lines.append("//")
    path.write_text("\n".join(lines) + "\n")
    return pol


def substring_aligner(reference_pol: str):
    """Fake MAFFT for sequences that are exact substrings of the reference."""

    def runner(ref_path: Path, frag_path: Path, threads: int) -> str:
        ((ref_name, ref_seq),) = read_fasta(ref_path)
        out = [f">{ref_name}\n{ref_seq}\n"]
        for name, seq in read_fasta(frag_path):
            offset = ref_seq.find(seq)
            assert offset >= 0, "test sequence must be an exact reference substring"
            row = "-" * offset + seq + "-" * (len(ref_seq) - offset - len(seq))
            out.append(f">{name}\n{row}\n")
        return "".join(out)

    return runner


@pytest.fixture
def workspace(config_path, tmp_path):
    pol = _hxb2_like_genbank(tmp_path / "hxb2.gb")
    settings = load_settings(
        config_path,
        "dev",
        {
            "paths": {
                "raw": str(tmp_path / "raw"),
                "interim": str(tmp_path / "interim"),
                "processed": str(tmp_path / "processed"),
                "reference": str(tmp_path),
                "results": str(tmp_path / "results"),
                "auspice": str(tmp_path / "auspice"),
                "reports": str(tmp_path / "reports"),
            },
            "align": {"chunk_size": 2},
        },
    )
    return settings, pol, tmp_path


def _seed_input(tmp_path, pol, count=5):
    strains = [(f"S{i}", pol[i * 30 : i * 30 + 900]) for i in range(count)]
    write_fasta(tmp_path / "interim" / "kenya" / "seqs.dedup.fasta", strains)
    return strains


def test_every_sequence_is_projected_in_order_across_chunks(workspace):
    settings, pol, tmp = workspace
    strains = _seed_input(tmp, pol, count=5)  # chunk_size 2 -> 3 chunks
    summary = run_align(settings, runner=substring_aligner(pol))

    assert summary["n_sequences"] == 5 and summary["n_chunks"] == 3
    rows = list(read_fasta(tmp / "interim" / "kenya" / "aln.pol.unmasked.fasta"))
    assert [name for name, _ in rows] == [name for name, _ in strains]
    assert {len(seq) for _, seq in rows} == {POL_END - POL_START + 1}
    # S1 starts 30 bases into pol and runs 900 bases
    s1 = dict(rows)["S1"]
    assert s1[:30] == "-" * 30 and s1[30:930] == strains[1][1] and set(s1[930:]) == {"-"}


def test_insertions_are_written_to_the_side_table(workspace):
    settings, pol, tmp = workspace
    _seed_input(tmp, pol, count=1)

    def runner(ref_path, frag_path, threads):
        ((ref_name, ref_seq),) = read_fasta(ref_path)
        ((name, seq),) = read_fasta(frag_path)
        cut = 600  # the sequence carries 2 extra bases after reference position 600
        ref_row = ref_seq[:cut] + "--" + ref_seq[cut:]
        seq_row = seq[:cut] + "GG" + seq[cut:] + "-" * (len(ref_seq) - len(seq))
        return f">{ref_name}\n{ref_row}\n>{name}\n{seq_row}\n"

    # sequence must be long enough to reach position 600
    write_fasta(tmp / "interim" / "kenya" / "seqs.dedup.fasta", [("S0", pol[:900])])
    summary = run_align(settings, runner=runner)
    assert summary["n_with_insertions"] == 1
    (row,) = read_tsv(tmp / "interim" / "kenya" / "insertions.tsv")
    assert row == {"strain": "S0", "after_ref_pos": "600", "length": "2"}


def test_a_sequence_missing_from_the_alignment_aborts_the_run(workspace):
    settings, pol, tmp = workspace
    _seed_input(tmp, pol, count=2)

    def lossy(ref_path, frag_path, threads):
        ((ref_name, ref_seq),) = read_fasta(ref_path)
        (first, seq), *_ = read_fasta(frag_path)
        return f">{ref_name}\n{ref_seq}\n>{first}\n{seq + '-' * (len(ref_seq) - len(seq))}\n"

    with pytest.raises(ContractError, match="missing"):
        run_align(settings, runner=lossy)


def test_missing_input_points_at_the_previous_stage(workspace):
    settings, pol, _ = workspace
    with pytest.raises(ContractError, match="normalise"):
        run_align(settings, runner=substring_aligner(pol))


# ------------------------------------------------------------------ mafft helpers
def test_chunking_covers_everything_once():
    items = [(f"s{i}", "ACGT") for i in range(7)]
    chunks = list(chunked(items, 3))
    assert [len(c) for c in chunks] == [3, 3, 1]
    assert [x for c in chunks for x in c] == items
    with pytest.raises(ContractError):
        list(chunked(items, 0))


def test_parse_aligned_fasta_joins_wrapped_lines_and_rejects_duplicates():
    assert parse_aligned_fasta(">a desc\nAC\nGT\n>b\nTT--\n") == {"a": "ACGT", "b": "TT--"}
    with pytest.raises(ContractError, match="duplicate"):
        parse_aligned_fasta(">a\nAC\n>a\nGT\n")


def test_mismatched_row_lengths_are_rejected():
    def bad(ref_path, frag_path, threads):
        return ">ref\nACGT\n>s1\nAC\n"

    with pytest.raises(ContractError, match="same length"):
        align_chunk(("ref", "ACGT"), [("s1", "AC")], runner=bad)
