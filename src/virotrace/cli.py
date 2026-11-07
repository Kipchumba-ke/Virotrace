"""Virotrace command line interface.

Every stage is independently runnable so it can be debugged, tested and
demonstrated in isolation. Snakemake wires stages together; it does not contain
logic.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Annotated

import typer

from . import __version__
from .errors import ConfigError, ContractError
from .provenance import tool_versions
from .settings import load_settings

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="HIV-1 drug-resistance and lineage surveillance for Kenya and Africa.",
)
drm_app = typer.Typer(no_args_is_help=True, help="Drug resistance layer.")
phylo_app = typer.Typer(no_args_is_help=True, help="Phylogenomics layer.")
export_app = typer.Typer(no_args_is_help=True, help="Output assembly.")
validate_app = typer.Typer(no_args_is_help=True, help="Validation without execution.")
app.add_typer(drm_app, name="drm")
app.add_typer(phylo_app, name="phylo")
app.add_typer(export_app, name="export")
app.add_typer(validate_app, name="validate")

ConfigOpt = Annotated[Path, typer.Option("--config", help="Base config YAML.")]
ProfileOpt = Annotated[str | None, typer.Option("--profile", help="dev | kenya | africa")]
SubtypeOpt = Annotated[
    str | None, typer.Option("--subtype", help="Build for one called subtype (e.g. A1) only.")
]


def _not_implemented(stage: str) -> None:
    typer.secho(
        f"stage '{stage}' is scaffolded but not implemented yet "
        f"(see docs/01-technical-specification.md section 8 for the phase it lands in)",
        fg=typer.colors.YELLOW,
        err=True,
    )
    raise typer.Exit(code=2)


@app.command()
def info() -> None:
    """Print versions of Virotrace and every external binary it shells out to."""
    payload = {"virotrace": __version__, "python": sys.version.split()[0], **tool_versions()}
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))


@validate_app.command("config")
def validate_config(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """Resolve and validate configuration without running anything.

    Exists so a bad config fails in a second rather than forty minutes into a tree.
    """
    try:
        settings = load_settings(config, profile)
    except (ConfigError, ContractError) as exc:
        typer.secho(f"invalid configuration: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1) from exc

    for warning in settings.warnings():
        typer.secho(f"warning: {warning}", fg=typer.colors.YELLOW, err=True)

    typer.secho(
        f"config valid  scope={settings.scope}  hash={settings.config_hash()}"
        f"  schema_version={settings.schema_version}",
        fg=typer.colors.GREEN,
    )


# --------------------------------------------------------------------- stage stubs
def _fail(exc: Exception) -> None:
    typer.secho(f"error: {exc}", fg=typer.colors.RED, err=True)
    raise typer.Exit(code=1) from exc


@app.command()
def acquire(
    config: ConfigOpt = Path("config/config.yaml"),
    profile: ProfileOpt = None,
    max_records: Annotated[
        int | None,
        typer.Option("--max-records", help="Cap the fetch for a small trial run."),
    ] = None,
) -> None:
    """L0: fetch sequences into an immutable snapshot with a manifest."""
    from .acquire.stage import run_acquire

    settings = load_settings(config, profile)
    try:
        snapshot = run_acquire(settings, max_records=max_records)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.secho(f"snapshot written: {snapshot}", fg=typer.colors.GREEN)


@app.command()
def normalise(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L1: dedup, scope filter and validate records; prints the P1 corpus counts."""
    from .curate.stage import run_normalise

    settings = load_settings(config, profile)
    try:
        summary = run_normalise(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def qc(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: QC rules on the aligned sequences, with the exclusion ledger and corpus recount."""
    from .curate.qc_stage import run_qc

    settings = load_settings(config, profile)
    try:
        summary = run_qc(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def align(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: align to HXB2 pol with MAFFT, projected onto reference coordinates."""
    from .align.stage import run_align

    settings = load_settings(config, profile)
    try:
        summary = run_align(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def mask(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: gap-mask the QC-passed alignment for tree inference (never for DRM calling)."""
    from .curate.mask import run_mask

    settings = load_settings(config, profile)
    try:
        summary = run_mask(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def subtype(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: bootstrap-supported subtype placement; depositor label compared, not trusted."""
    from .curate.subtype import run_subtype

    settings = load_settings(config, profile)
    try:
        summary = run_subtype(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def recombination(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L2: breakpoint screen. Confirmed recombinants leave the tree, stay in DRM."""
    from .curate.recombination import run_recombination

    settings = load_settings(config, profile)
    try:
        summary = run_recombination(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command()
def subsample(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: deterministic stratified subsampling."""
    from .curate.subsample_stage import run_subsample

    settings = load_settings(config, profile)
    try:
        summary = run_subsample(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@drm_app.command("call")
def drm_call(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L3a: per-gene coverage and mutations at SDRM positions, from the pol alignment."""
    from .drm.stage import run_drm_call

    settings = load_settings(config, profile)
    try:
        summary = run_drm_call(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@drm_app.command("prevalence")
def drm_prevalence(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3a: per-class denominators, Wilson intervals, suppression, trend."""
    from .drm.prevalence_stage import run_prevalence

    settings = load_settings(config, profile)
    try:
        summary = run_prevalence(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@phylo_app.command("tree")
def phylo_tree(
    config: ConfigOpt = Path("config/config.yaml"),
    profile: ProfileOpt = None,
    subtype: SubtypeOpt = None,
) -> None:
    """L3b: maximum-likelihood tree."""
    from .phylo.stage import run_tree

    settings = load_settings(config, profile)
    try:
        summary = run_tree(settings, subtype=subtype)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@phylo_app.command("refine")
def phylo_refine(
    config: ConfigOpt = Path("config/config.yaml"),
    profile: ProfileOpt = None,
    subtype: SubtypeOpt = None,
) -> None:
    """L3b: time-resolved refinement behind the temporal-signal gate."""
    from .phylo.stage import run_refine

    settings = load_settings(config, profile)
    try:
        summary = run_refine(settings, subtype=subtype)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@phylo_app.command("subtrees")
def phylo_subtrees(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3b: a tree (and, if its own gate passes, a time tree) for each configured subtype."""
    from .phylo.stage import run_subtype_trees

    settings = load_settings(config, profile)
    try:
        summary = run_subtype_trees(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@phylo_app.command("traits")
def phylo_traits(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3b: ancestral trait reconstruction over country and region."""
    load_settings(config, profile)
    _not_implemented("phylo traits")


@phylo_app.command("clusters")
def phylo_clusters(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3b: TN93 clustering as a descriptive population statistic."""
    from .phylo.clusters import run_clusters

    settings = load_settings(config, profile)
    try:
        summary = run_clusters(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@export_app.command("auspice")
def export_auspice(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L4: Auspice v2 JSON with DRM colourings."""
    from .export.auspice import run_export

    settings = load_settings(config, profile)
    try:
        summary = run_export(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@export_app.command("report")
def export_report(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L4: the surveillance report, limitations section included by construction."""
    from .export.report import run_report

    settings = load_settings(config, profile)
    try:
        summary = run_report(settings)
    except (ConfigError, ContractError) as exc:
        _fail(exc)
        return
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


def main() -> None:
    app()


if __name__ == "__main__":
    main()
