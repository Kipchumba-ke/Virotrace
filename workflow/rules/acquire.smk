"""L0 -- acquisition into an immutable snapshot.

data/raw is never declared as an output by any other rule; a CI check enforces that,
which is what makes "immutable snapshot" true rather than aspirational.
"""


rule acquire:
    output:
        seqs=f"{PATHS['raw']}/{SCOPE}/sequences.fasta",
        meta=f"{PATHS['raw']}/{SCOPE}/records.tsv",
        manifest=f"{PATHS['raw']}/{SCOPE}/manifest.json",
    log:
        f"logs/{SCOPE}/acquire.log",
    threads: 1
    shell:
        "virotrace acquire {CONFIG_ARGS} > {log} 2>&1"
