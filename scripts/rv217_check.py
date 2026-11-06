"""Do the RV217 individuals called C look like C or like the A1 individuals?

For one sequence per RV217 individual (the one with most comparable sites), print the
call, the distance to the nearest pure-panel reference of each label, and the nearest
*other* individual with that individual's call. Run from the repo root:
    python scripts/rv217_check.py kenya
"""

import re
import sys

import numpy as np

from virotrace.curate.subtype import encode, is_pure, panel_label
from virotrace.io.fasta import read_fasta
from virotrace.io.tables import read_tsv

scope = sys.argv[1] if len(sys.argv) > 1 else "kenya"
d = f"data/interim/{scope}"
pat = re.compile(r"^(RV217_\d+)_")
calls = {r["strain"]: r["subtype_called"] for r in read_tsv(f"{d}/subtype_calls.tsv")}
iso = {r["strain"]: r.get("isolate", "") for r in read_tsv(f"{d}/metadata.valid.tsv")}

best = {}
for name, seq in read_fasta(f"{d}/aln.pol.unmasked.fasta"):
    m = pat.match(iso.get(name, ""))
    if not m or name not in calls:
        continue
    e = encode(seq)
    n = int((e < 4).sum())
    if m.group(1) not in best or n > best[m.group(1)][2]:
        best[m.group(1)] = (name, e, n)

panel = [(h, encode(s)) for h, s in read_fasta(f"{d}/panel.pol.fasta") if is_pure(panel_label(h))]


def dist(a, b):
    ok = (a < 4) & (b < 4)
    return float((a[ok] != b[ok]).mean()) if ok.sum() >= 300 else float("nan")


ids = sorted(best)
labels = sorted({panel_label(h) for h, _ in panel})
print("individual", "call", *labels, "nearest_other(call)", sep="\t")
for i in ids:
    name, e, _ = best[i]
    per = []
    for lab in labels:
        ds = [dist(e, p) for h, p in panel if panel_label(h) == lab]
        per.append(np.nanmin(ds))
    others = sorted((dist(e, best[j][1]), j) for j in ids if j != i)
    dj, j = others[0]
    print(
        i, calls[name], *[f"{x:.3f}" for x in per], f"{j} {dj:.3f} ({calls[best[j][0]]})", sep="\t"
    )
