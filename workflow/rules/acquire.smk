"""L0 -- acquisition into an immutable snapshot.

Snapshot directories are named dynamically (`<date>_<query hash>`), which Snakemake
cannot declare as a fixed output. So the rule's declared output is a small pointer
file, `data/raw/<scope>.current`, holding the snapshot id. Downstream rules depend on
the pointer; the normalise stage resolves it to the snapshot directory (or uses
`snapshot_id` from the config when pinned, for a reproducible rerun).

data/raw is never written by any rule other than this one; a CI check enforces that.
"""


rule acquire:
    output:
        pointer=f"{PATHS['raw']}/{SCOPE}.current",
    log:
        f"logs/{SCOPE}/acquire.log",
    threads: 1
    shell:
        "virotrace acquire {CONFIG_ARGS} > {log} 2>&1"
