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


SUBTYPE_TREES = config["phylo"]["refine"]["subtype_trees"]


if SUBTYPE_TREES:  # dev has none; an empty rule would have no outputs

    rule subtrees:
        """One tree per configured subtype, each with its own temporal gate.

        A pooled multi-subtype tree rarely has a usable clock; a subtype's own tree often
        does. Each subtype gets IQ-TREE, the gate, and a time tree only if its gate passes.
        """
        input:
            tree_aln=rules.subsample.output.tree_aln,
            meta=rules.subsample.output.meta,
        output:
            tree=expand(f"{RESULTS}/phylo/{{subtype}}/tree.refined.nwk", subtype=SUBTYPE_TREES),
            node_data=expand(f"{RESULTS}/phylo/{{subtype}}/branch_lengths.json", subtype=SUBTYPE_TREES),
            gate=expand(f"{RESULTS}/phylo/{{subtype}}/temporal_gate.json", subtype=SUBTYPE_TREES),
        log:
            f"logs/{SCOPE}/subtrees.log",
        threads: config["threads"]
        shell:
            "virotrace phylo subtrees {CONFIG_ARGS} > {log} 2>&1"


if config["phylo"]["traits"]["enabled"]:

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
        summary=f"{RESULTS}/phylo/clusters_summary.json",
        near_identical=f"{RESULTS}/qc/near_identical_pairs.tsv",
    log:
        f"logs/{SCOPE}/clusters.log",
    threads: config["threads"]
    shell:
        "virotrace phylo clusters {CONFIG_ARGS} > {log} 2>&1"
