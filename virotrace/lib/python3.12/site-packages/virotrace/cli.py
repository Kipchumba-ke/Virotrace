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
@app.command()
def acquire(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L0: fetch sequences into an immutable snapshot with a manifest."""
    load_settings(config, profile)
    _not_implemented("acquire")


@app.command()
def normalise(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L1: source records -> the validated metadata contract."""
    load_settings(config, profile)
    _not_implemented("normalise")


@app.command()
def qc(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: QC-1..QC-7 with the exclusion ledger."""
    load_settings(config, profile)
    _not_implemented("qc")


@app.command()
def align(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: reference-anchored, codon-aware alignment to HXB2."""
    load_settings(config, profile)
    _not_implemented("align")


@app.command()
def subtype(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: two independent subtype calls, disagreement surfaced not resolved."""
    load_settings(config, profile)
    _not_implemented("subtype")


@app.command()
def recombination(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L2: breakpoint screen. Confirmed recombinants leave the tree, stay in DRM."""
    load_settings(config, profile)
    _not_implemented("recombination")


@app.command()
def subsample(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L2: deterministic stratified subsampling."""
    load_settings(config, profile)
    _not_implemented("subsample")


@drm_app.command("call")
def drm_call(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L3a: mutation extraction and HIVdb interpretation."""
    load_settings(config, profile)
    _not_implemented("drm call")


@drm_app.command("prevalence")
def drm_prevalence(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3a: per-class denominators, Wilson intervals, suppression, trend."""
    load_settings(config, profile)
    _not_implemented("drm prevalence")


@phylo_app.command("tree")
def phylo_tree(config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None) -> None:
    """L3b: maximum-likelihood tree."""
    load_settings(config, profile)
    _not_implemented("phylo tree")


@phylo_app.command("refine")
def phylo_refine(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L3b: time-resolved refinement behind the temporal-signal gate."""
    load_settings(config, profile)
    _not_implemented("phylo refine")


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
    load_settings(config, profile)
    _not_implemented("phylo clusters")


@export_app.command("auspice")
def export_auspice(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L4: Auspice v2 JSON with DRM colourings."""
    load_settings(config, profile)
    _not_implemented("export auspice")


@export_app.command("report")
def export_report(
    config: ConfigOpt = Path("config/config.yaml"), profile: ProfileOpt = None
) -> None:
    """L4: the surveillance report, limitations section included by construction."""
    load_settings(config, profile)
    _not_implemented("export report")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
