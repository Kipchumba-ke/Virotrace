#!/usr/bin/env python
"""Build config/reference/sdrm_who.tsv from a transcribed copy of the WHO table.

    python scripts/sdrm_from_text.py sdrm_table.txt --version 2009 \
        --source "Bennett et al. 2009, PLoS ONE 4(3):e4724"

Prints counts per class and any wild-type residue that disagrees with HXB2.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

from virotrace.drm.reference import gene_proteins, load_reference_sequence
from virotrace.drm.sdrm_import import check_wild_types, parse_sdrm_text, render_tsv
from virotrace.errors import ContractError
from virotrace.settings import load_settings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("text", type=Path, help="transcribed table, one class per line")
    parser.add_argument("--version", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--out", type=Path, default=Path("config/reference/sdrm_who.tsv"))
    args = parser.parse_args()

    text = args.text.read_text()
    try:
        rows = parse_sdrm_text(text)
    except ContractError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    counts = collections.Counter(r.drug_class for r in rows)
    print("positions per class:", dict(counts))

    settings = load_settings(args.config)
    record = args.config.parent / "reference" / settings.reference.record
    if record.exists():
        genes = settings.reference.genes
        proteins = gene_proteins(
            load_reference_sequence(record),
            {"PR": tuple(genes.pr), "RT": tuple(genes.rt), "IN": tuple(genes.in_)},
        )
        problems = check_wild_types(rows, proteins)
        if problems:
            print("WARNING: wild-type residues that differ from HXB2 (re-read the source):")
            for problem in problems:
                print("  ", problem)
        else:
            print("all wild-type residues agree with HXB2")
    else:
        print(f"note: {record} not found, wild-type check skipped")

    args.out.write_text(render_tsv(rows, source=args.source, version=args.version, input_text=text))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
