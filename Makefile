.DEFAULT_GOAL := help
SHELL := /bin/bash
PROFILE ?= dev
CORES ?= 4

.PHONY: help setup lint type test test-all smoke dag dev kenya africa refresh-golden clean clean-all

help: ## Show this help
	@grep -hE '^[a-z-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

setup: ## Create the conda environment and install the package
	mamba env create -f environment.yml || conda env create -f environment.yml
	@echo "now run: conda activate virotrace"

lint: ## Lint and format check
	ruff check src tests
	ruff format --check src tests

type: ## Static type check
	mypy

test: ## Fast suite (unit, contract, integration) -- target under 60s
	pytest tests/unit tests/contract tests/integration

test-all: ## Everything including the smoke test
	pytest

smoke: ## Full dev-profile DAG end to end
	pytest -m smoke

dag: ## Validate the DAG for every profile without executing anything
	snakemake -s workflow/Snakefile --configfile config/config.yaml --dry-run --cores 1
	virotrace validate config --config config/config.yaml --profile dev
	virotrace validate config --config config/config.yaml --profile kenya
	virotrace validate config --config config/config.yaml --profile africa

dev: ## Run the toy build (minutes, no network)
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 2 \
		--config scope=kenya

kenya: ## Run the primary Kenya build
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores $(CORES) \
		--config scope=kenya

africa: ## Run the Africa context build (long; resumable)
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores $(CORES) \
		--config scope=africa --rerun-incomplete

refresh-golden: ## Regenerate regression baselines -- the diff must be explained in the PR
	pytest tests/regression --snapshot-update || true
	@echo "review the diff before committing; an unexplained refresh blesses a regression"

clean: ## Remove derived results, keep raw snapshots
	rm -rf data/interim data/processed results auspice reports logs .snakemake

clean-all: clean ## Also remove raw snapshots
	rm -rf data/raw
