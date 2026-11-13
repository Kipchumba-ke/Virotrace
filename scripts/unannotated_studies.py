"""List the largest studies that still have no treatment annotation.

Run from the repo root after `virotrace qc` has run for the scope:
    python scripts/unannotated_studies.py africa [top_n]

Prints, per study: QC-passed sequences, countries, years, the study title and the keyword
hint distribution. The hint is for orientation only; the annotation must come from the paper.
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    scope = sys.argv[1] if len(sys.argv) > 1 else "africa"
    top = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    interim = Path("data/interim") / scope
    passed = {r["strain"] for r in rows(interim / "qc.passed.tsv")}
    meta = [r for r in rows(interim / "metadata.valid.tsv") if r["strain"] in passed]
    ann = Path("config/reference/study_annotations.tsv")
    done = set()
    if ann.exists():
        with ann.open(newline="") as handle:
            lines = (ln for ln in handle if not ln.lstrip().startswith("#"))
            done = {r["study_key"].strip() for r in csv.DictReader(lines, delimiter="\t")}
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in meta:
        groups[r.get("study_key", "") or "unknown_study"].append(r)
    total = len(meta)
    todo = sorted(((k, v) for k, v in groups.items() if k not in done), key=lambda kv: -len(kv[1]))
    covered = sum(len(v) for k, v in groups.items() if k in done)
    print(
        f"{total} sequences; {covered} annotated ({covered / total:.1%}); {len(todo)} studies open"
    )
    running = covered
    for key, members in todo[:top]:
        running += len(members)
        countries = Counter(m.get("country", "") for m in members).most_common(3)
        years = sorted(int(m["year"]) for m in members if m.get("year", "").strip().isdigit())
        span = f"{years[0]}-{years[-1]}" if years else "no years"
        hints = dict(Counter(m.get("treatment_status", "") or "-" for m in members))
        title = next((m.get("study_title", "") for m in members if m.get("study_title")), "")
        print(f"\n{key}\n  n={len(members)} (cumulative {running / total:.1%})  {span}")
        print(f"  countries: {countries}\n  hint: {hints}\n  title: {title[:140]}")


if __name__ == "__main__":
    main()
