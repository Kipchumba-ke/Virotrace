"""Print full authors, journal/submitter line, accessions and notes for chosen study keys.

Reads the cached GenBank flat files (data/raw/.fetch_cache/*.gb), which still carry the
REFERENCE block that the metadata drops. Use it to identify "Direct Submission" studies:

    python scripts/study_authors.py "REF:Ulenga|direct-submission" "REF:Fokam|direct-submission"
    python scripts/study_authors.py --prefix REF: --only-direct   # every direct submission
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

from virotrace.acquire.genbank_parse import parse_genbank
from virotrace.curate.study import study_key


def main() -> None:
    args = sys.argv[1:]
    direct = "--only-direct" in args
    wanted = {a for a in args if not a.startswith("--")}
    cache = Path("data/raw/.fetch_cache")
    groups: dict[str, list] = defaultdict(list)
    for gb in sorted(cache.glob("*.gb")):
        for rec in parse_genbank(gb.read_text(), source="genbank"):
            key = study_key(rec.fields)
            if key in wanted or (direct and key.endswith("|direct-submission")):
                groups[key].append(rec)
    for key, recs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        print(f"\n{key}  n={len(recs)}")
        for field in ("study_authors", "study_journal"):
            for value, count in Counter(r.fields.get(field, "") for r in recs).most_common(4):
                print(f"  {field} x{count}: {value[:300]}")
        for field in ("isolation_source", "note", "host", "definition"):
            vals = Counter(r.fields.get(field, "") for r in recs if r.fields.get(field, ""))
            if vals:
                print(f"  {field}: {vals.most_common(3)}")
        print("  accessions:", ", ".join(r.accession for r in recs[:6]))


if __name__ == "__main__":
    main()
