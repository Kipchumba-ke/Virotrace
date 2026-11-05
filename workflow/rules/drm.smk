"""L3a -- drug resistance layer.

Consumes the unmasked alignment: a masked position must become a no-call, never a
wild-type call. This layer has no dependency on the phylo layer.
"""


rule drm_prevalence:
    # Counts come from the QC-passed, clone-policy-filtered metadata, not from the
    # subsampled tree set: subsampling is for tree inference and must not shape prevalence.
    input:
        table=rules.drm_call.output.table,
        coverage=rules.drm_call.output.coverage,
        meta=rules.normalise.output.valid,
        passed=rules.qc.output.passed,
        subtypes=rules.subtype.output.calls,
    output:
        prevalence=f"{RESULTS}/drm/prevalence.tsv",
        trend=f"{RESULTS}/drm/trend_tests.tsv",
        summary=f"{RESULTS}/drm/prevalence_summary.json",
        studies=f"{RESULTS}/drm/study_summary.tsv",
    log:
        f"logs/{SCOPE}/drm_prevalence.log",
    shell:
        "virotrace drm prevalence {CONFIG_ARGS} > {log} 2>&1"
