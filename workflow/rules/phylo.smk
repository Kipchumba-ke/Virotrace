"""L3b -- phylogenomics layer.

Consumes the masked tree alignment and the recombination flags. Joins the DRM layer
only through drm_table.tsv, which is a file contract, not a code dependency.
"""


rule tree:
    input:
        aln=rules.subsample.output.tree_aln,
    output:
        tree=f"{RESULTS}/phylo/tree.nwk",
        iqtree_log=f"{RESULTS}/phylo/iqtree.log",
    log:
        f"logs/{SCOPE}/tree.log",
    threads: config["threads"]
    shell:
        "virotrace phylo tree {CONFIG_ARGS} > {log} 2>&1"


rule refine:
    """Time-resolved refinement behind the temporal-signal gate.

    If root-to-tip regression fails the configured threshold, this rule still
    produces its outputs -- carrying the divergence tree and a recorded reason
    rather than a time tree nobody should believe.
    """
    input:
        tree=rules.tree.output.tree,
        aln=rules.subsample.output.tree_aln,
        meta=rules.subsample.output.meta,
    output:
        tree=f"{RESULTS}/phylo/tree.refined.nwk",
        node_data=f"{RESULTS}/phylo/branch_lengths.json",
        gate=f"{RESULTS}/phylo/temporal_gate.json",
    log:
        f"logs/{SCOPE}/refine.log",
    threads: config["threads"]
    shell:
        "virotrace phylo refine {CONFIG_ARGS} > {log} 2>&1"


rule traits:
    input:
        tree=rules.refine.output.tree,
        meta=rules.subsample.output.meta,
    output:
        traits=f"{RESULTS}/phylo/traits.json",
    log:
        f"logs/{SCOPE}/traits.log",
    shell:
        "virotrace phylo traits {CONFIG_ARGS} > {log} 2>&1"


rule clusters:
    input:
        aln=rules.subsample.output.tree_aln,
        meta=rules.subsample.output.meta,
        drm=rules.drm_call.output.table,
    output:
        clusters=f"{RESULTS}/phylo/clusters.tsv",
    log:
        f"logs/{SCOPE}/clusters.log",
    threads: config["threads"]
    shell:
        "virotrace phylo clusters {CONFIG_ARGS} > {log} 2>&1"
