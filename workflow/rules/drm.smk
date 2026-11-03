"""L3a -- drug resistance layer.

Consumes the unmasked alignment: a masked position must become a no-call, never a
wild-type call. This layer has no dependency on the phylo layer.
"""


rule drm_call:
    input:
        aln=rules.subsample.output.drm_aln,
        meta=rules.subsample.output.meta,
        sdrm=f"{PATHS['reference']}/{config['drm']['sdrm']['list']}",
    output:
        table=f"{RESULTS}/drm/drm_table.tsv",
        coverage=f"{RESULTS}/drm/gene_coverage.tsv",
    log:
        f"logs/{SCOPE}/drm_call.log",
    threads: config["threads"]
    shell:
        "virotrace drm call {CONFIG_ARGS} > {log} 2>&1"


rule drm_prevalence:
    input:
        table=rules.drm_call.output.table,
        coverage=rules.drm_call.output.coverage,
        meta=rules.subsample.output.meta,
    output:
        prevalence=f"{RESULTS}/drm/prevalence.tsv",
        trend=f"{RESULTS}/drm/trend_tests.tsv",
    log:
        f"logs/{SCOPE}/drm_prevalence.log",
    shell:
        "virotrace drm prevalence {CONFIG_ARGS} > {log} 2>&1"
