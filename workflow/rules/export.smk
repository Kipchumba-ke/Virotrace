"""L4 -- export. Assembles outputs for humans; computes nothing new."""


def auspice_inputs(wildcards):
    required = {
        "meta": rules.subsample.output.meta,
        "drm": rules.drm_call.output.table,
    }
    if PHYLO_ENABLED:
        required.update(
            {
                "tree": rules.refine.output.tree,
                "node_data": rules.refine.output.node_data,
                "traits": rules.traits.output.traits,
                "gate": rules.refine.output.gate,
            }
        )
    if PHYLO_ENABLED and config["phylo"]["clusters"]["enabled"]:
        required["clusters"] = rules.clusters.output.clusters
    return required


rule resolved_config:
    """Write the fully resolved config beside the results.

    The config in the repo is the input; this one is the record of what actually ran.
    """
    output:
        resolved=f"{RESULTS}/resolved_config.yaml",
    shell:
        "virotrace validate config {CONFIG_ARGS} && "
        "mkdir -p $(dirname {output.resolved}) && "
        "cp config/config.yaml {output.resolved}"


rule export_auspice:
    input:
        unpack(auspice_inputs),
    output:
        auspice=f"{PATHS['auspice']}/virotrace_{SCOPE}_pol.json",
    log:
        f"logs/{SCOPE}/export_auspice.log",
    shell:
        "virotrace export auspice {CONFIG_ARGS} > {log} 2>&1"


rule export_report:
    input:
        prevalence=rules.drm_prevalence.output.prevalence,
        trend=rules.drm_prevalence.output.trend,
        ledger=rules.qc.output.ledger,
        strata=rules.subsample.output.strata,
        resolved=rules.resolved_config.output.resolved,
    output:
        report=f"{PATHS['reports']}/{SCOPE}/index.html",
    log:
        f"logs/{SCOPE}/export_report.log",
    shell:
        "virotrace export report {CONFIG_ARGS} > {log} 2>&1"
