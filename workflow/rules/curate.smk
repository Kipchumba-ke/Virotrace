"""L1 + L2 -- normalisation, QC, alignment, subtyping, recombination, subsampling.

The metadata.valid.tsv / aln.final.fasta pair is the contract both analysis layers
consume. Nothing downstream of here reaches back past it.
"""


rule normalise:
    input:
        pointer=rules.acquire.output.pointer,
    output:
        valid=f"{INTERIM}/metadata.valid.tsv",
        seqs=f"{INTERIM}/seqs.dedup.fasta",
        quarantine=f"{RESULTS}/qc/quarantine.tsv",
        dedup_log=f"{RESULTS}/qc/dedup_log.tsv",
        out_of_scope=f"{RESULTS}/qc/out_of_scope.tsv",
        counts=f"{RESULTS}/qc/corpus_counts.json",
    log:
        f"logs/{SCOPE}/normalise.log",
    shell:
        "virotrace normalise {CONFIG_ARGS} > {log} 2>&1"


rule align:
    # Alignment comes BEFORE QC: coverage, internal stops and frameshifts (QC-1, 3, 4)
    # are defined on HXB2 coordinates, which only exist once a sequence is aligned.
    input:
        seqs=rules.normalise.output.seqs,
        reference=f"{PATHS['reference']}/{config['reference']['record']}",
    output:
        # Rectangular alignment in HXB2 pol coordinates. Never masked: masking is an
        # acceptable tree convenience and an unacceptable resistance-calling shortcut.
        drm_aln=f"{INTERIM}/aln.pol.unmasked.fasta",
        insertions=f"{INTERIM}/insertions.tsv",
        summary=f"{RESULTS}/qc/align_summary.json",
    log:
        f"logs/{SCOPE}/align.log",
    threads: config["threads"]
    shell:
        "virotrace align {CONFIG_ARGS} > {log} 2>&1"


rule drm_call:
    # Lives beside align (not in drm.smk) because QC consumes its coverage and identity
    # table, and Snakemake needs a rule defined before it is referenced.
    # Reads the aligned pol sequences directly. It does not wait for subtyping,
    # recombination or subsampling: the resistance layer is independent of the phylo
    # layer, and which rows to count is decided later against the coverage table.
    input:
        aln=rules.align.output.drm_aln,
        insertions=rules.align.output.insertions,
        sdrm=f"{PATHS['reference']}/{config['drm']['sdrm']['list']}",
        reference=f"{PATHS['reference']}/{config['reference']['record']}",
    output:
        table=f"{RESULTS}/drm/drm_table.tsv",
        coverage=f"{RESULTS}/drm/gene_coverage.tsv",
        summary=f"{RESULTS}/drm/drm_summary.json",
    log:
        f"logs/{SCOPE}/drm_call.log",
    threads: config["threads"]
    shell:
        "virotrace drm call {CONFIG_ARGS} > {log} 2>&1"


rule qc:
    input:
        aln=rules.align.output.drm_aln,
        coverage=rules.drm_call.output.coverage,
        meta=rules.normalise.output.valid,
        lab_strains=f"{PATHS['reference']}/{config['qc']['lab_strain_list']}",
    output:
        passed=f"{INTERIM}/qc.passed.tsv",
        ledger=f"{RESULTS}/qc/excluded.tsv",
        summary=f"{RESULTS}/qc/qc_summary.json",
    log:
        f"logs/{SCOPE}/qc.log",
    shell:
        "virotrace qc {CONFIG_ARGS} > {log} 2>&1"


rule mask_alignment:
    input:
        aln=rules.align.output.drm_aln,
        passed=rules.qc.output.passed,
    output:
        tree_aln=f"{INTERIM}/aln.pol.masked.fasta",
    log:
        f"logs/{SCOPE}/mask_alignment.log",
    shell:
        "virotrace mask {CONFIG_ARGS} > {log} 2>&1"


rule subtype:
    # Reads the unmasked alignment: every QC-passed sequence is subtyped, including the
    # ones that will not reach the tree, because the resistance layer reports by subtype.
    input:
        aln=rules.align.output.drm_aln,
        passed=rules.qc.output.passed,
        meta=rules.normalise.output.valid,
        panel=f"{PATHS['reference']}/{config['subtype']['reference_panel']}",
    output:
        calls=f"{INTERIM}/subtype_calls.tsv",
        panel_aln=f"{INTERIM}/panel.pol.fasta",
    log:
        f"logs/{SCOPE}/subtype.log",
    threads: config["threads"]
    shell:
        "virotrace subtype {CONFIG_ARGS} > {log} 2>&1"


rule recombination:
    input:
        aln=rules.align.output.drm_aln,
        panel_aln=rules.subtype.output.panel_aln,
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
        tree_aln=rules.mask_alignment.output.tree_aln,
    output:
        tree_aln=f"{PROCESSED}/aln.tree.fasta",
        meta=f"{PROCESSED}/metadata.final.tsv",
        strata=f"{RESULTS}/qc/subsample_strata.tsv",
    log:
        f"logs/{SCOPE}/subsample.log",
    shell:
        "virotrace subsample {CONFIG_ARGS} > {log} 2>&1"
