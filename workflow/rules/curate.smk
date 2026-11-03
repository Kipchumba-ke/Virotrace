"""L1 + L2 -- normalisation, QC, alignment, subtyping, recombination, subsampling.

The metadata.valid.tsv / aln.final.fasta pair is the contract both analysis layers
consume. Nothing downstream of here reaches back past it.
"""


rule normalise:
    input:
        meta=rules.acquire.output.meta,
        manifest=rules.acquire.output.manifest,
    output:
        valid=f"{INTERIM}/metadata.valid.tsv",
        quarantine=f"{RESULTS}/qc/quarantine.tsv",
    log:
        f"logs/{SCOPE}/normalise.log",
    shell:
        "virotrace normalise {CONFIG_ARGS} > {log} 2>&1"


rule qc:
    input:
        seqs=rules.acquire.output.seqs,
        meta=rules.normalise.output.valid,
    output:
        seqs=f"{INTERIM}/seqs.qc.fasta",
        ledger=f"{RESULTS}/qc/excluded.tsv",
    log:
        f"logs/{SCOPE}/qc.log",
    shell:
        "virotrace qc {CONFIG_ARGS} > {log} 2>&1"


rule align:
    input:
        seqs=rules.qc.output.seqs,
        reference=f"{PATHS['reference']}/{config['reference']['record']}",
    output:
        # Tree alignment (gap-masked) and DRM alignment (never masked) are kept
        # separate on purpose: masking is an acceptable tree convenience and an
        # unacceptable resistance-calling shortcut.
        tree_aln=f"{INTERIM}/aln.pol.masked.fasta",
        drm_aln=f"{INTERIM}/aln.pol.unmasked.fasta",
    log:
        f"logs/{SCOPE}/align.log",
    threads: config["threads"]
    shell:
        "virotrace align {CONFIG_ARGS} > {log} 2>&1"


rule subtype:
    input:
        aln=rules.align.output.tree_aln,
        panel=f"{PATHS['reference']}/{config['subtype']['reference_panel']}",
    output:
        calls=f"{INTERIM}/subtype_calls.tsv",
    log:
        f"logs/{SCOPE}/subtype.log",
    threads: config["threads"]
    shell:
        "virotrace subtype {CONFIG_ARGS} > {log} 2>&1"


rule recombination:
    input:
        aln=rules.align.output.tree_aln,
        subtypes=rules.subtype.output.calls,
    output:
        flags=f"{INTERIM}/recombination.tsv",
    log:
        f"logs/{SCOPE}/recombination.log",
    threads: config["threads"]
    shell:
        "virotrace recombination {CONFIG_ARGS} > {log} 2>&1"


rule subsample:
    input:
        meta=rules.normalise.output.valid,
        subtypes=rules.subtype.output.calls,
        recomb=rules.recombination.output.flags,
        tree_aln=rules.align.output.tree_aln,
        drm_aln=rules.align.output.drm_aln,
    output:
        tree_aln=f"{PROCESSED}/aln.tree.fasta",
        drm_aln=f"{PROCESSED}/aln.drm.fasta",
        meta=f"{PROCESSED}/metadata.final.tsv",
        strata=f"{RESULTS}/qc/subsample_strata.tsv",
    log:
        f"logs/{SCOPE}/subsample.log",
    shell:
        "virotrace subsample {CONFIG_ARGS} > {log} 2>&1"
