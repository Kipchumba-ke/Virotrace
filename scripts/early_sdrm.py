"""Show QC-passed individuals sampled in or before a cut-off year that carry an SDRM.

Pre-ART samples (before about 1996-2000) carrying NRTI SDRMs are either mislabelled dates,
prior exposure, or mutations that are polymorphic rather than selected. Run from the repo root:
    python scripts/early_sdrm.py africa [max_year=1989] [drug_class=NRTI]
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from virotrace.drm.classes import drug_class_of


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    scope = sys.argv[1] if len(sys.argv) > 1 else "africa"
    max_year = int(sys.argv[2]) if len(sys.argv) > 2 else 1989
    wanted = sys.argv[3] if len(sys.argv) > 3 else "NRTI"
    passed = {r["strain"] for r in rows(Path("data/interim") / scope / "qc.passed.tsv")}
    meta = {
        r["strain"]: r
        for r in rows(Path("data/interim") / scope / "metadata.valid.tsv")
        if r["strain"] in passed
    }
    early = {
        s
        for s, r in meta.items()
        if r.get("year", "").strip().isdigit() and int(r["year"]) <= max_year
    }
    hits: dict[str, list[str]] = defaultdict(list)
    for r in rows(Path("results") / scope / "drm" / "drm_table.tsv"):
        if (
            r["strain"] in early
            and r["is_sdrm"].strip().lower() == "true"
            and drug_class_of(r["gene"], int(r["position"])) == wanted
        ):
            hits[r["strain"]].append(f"{r['gene']}{r['position']}")
    print(f"{len(early)} QC-passed individuals dated <= {max_year}; {len(hits)} with {wanted} SDRM")
    for strain, muts in sorted(hits.items()):
        m = meta[strain]
        print(f"\n{strain}  {m.get('country', '')} {m.get('year', '')}  {sorted(muts)}")
        print(f"  study: {m.get('study_key', '')}  | {m.get('study_title', '')[:110]}")
        print(f"  date_raw/precision: {m.get('collection_date', '')} {m.get('date_precision', '')}")


if __name__ == "__main__":
    main()
