#!/usr/bin/env bash
# One-time setup: build the conda env and install DoubletFinder from GitHub.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME=scrna-demo

# Pinned DoubletFinder commit — resolved from
# https://github.com/chris-mcginnis-ucsf/DoubletFinder HEAD.
# Bump intentionally; do not rely on implicit 'main'.
DOUBLETFINDER_SHA=1b244d8f0d54b4b1cb4365639931bbb16f01e1cd

# Mirror all setup output to setuplog.txt alongside the script.
SETUP_LOG="$HERE/setuplog.txt"
: > "$SETUP_LOG"
exec > >(tee -a "$SETUP_LOG") 2>&1

echo "[setup.sh] Started $(date '+%Y-%m-%d %H:%M:%S')"
echo "[setup.sh] Log: $SETUP_LOG"

# Prefer mamba for speed; fall back to conda.
if command -v mamba >/dev/null 2>&1; then
    SOLVER=mamba
else
    SOLVER=conda
fi

# Bioconductor packages (SingleR, celldex, etc.) have incomplete native
# osx-arm64 coverage; force osx-64 on Darwin when conda is configured for
# arm64 or when this shell is already running under x86_64/Rosetta.
CONDA_PLATFORM="$(conda info --json 2>/dev/null | python -c 'import json,sys; print(json.load(sys.stdin).get("platform", ""))' 2>/dev/null || true)"
if [ "$(uname -s)" = "Darwin" ] && { [ "$(uname -m)" = "x86_64" ] || [ "$CONDA_PLATFORM" = "osx-arm64" ]; }; then
    export CONDA_SUBDIR=osx-64
    echo "==> Darwin detected with Bioconductor constraints; creating osx-64 env"
fi

echo "==> Creating env '$ENV_NAME' with $SOLVER (may take 5-10 min)"
$SOLVER env create -f "$HERE/environment.yml" || {
    echo "Env already exists? Attempting update instead..."
    $SOLVER env update -f "$HERE/environment.yml"
}

# Pin the subdir inside the env so future installs stay on osx-64.
if [ "${CONDA_SUBDIR:-}" = "osx-64" ]; then
    conda config --env --set subdir osx-64 --file \
        "$(conda info --base)/envs/$ENV_NAME/.condarc" 2>/dev/null || \
        conda run -n "$ENV_NAME" conda config --env --set subdir osx-64
fi

# Use the env's Rscript by absolute path; `conda run` can fall back to a
# system R on PATH and install packages into the wrong library.
ENV_PREFIX="$(conda info --base)/envs/$ENV_NAME"
ENV_RSCRIPT="$ENV_PREFIX/bin/Rscript"

echo "==> Installing DoubletFinder from GitHub (pinned SHA: $DOUBLETFINDER_SHA)"
"$ENV_RSCRIPT" -e \
    "if (!requireNamespace('DoubletFinder', quietly=TRUE))
         remotes::install_github('chris-mcginnis-ucsf/DoubletFinder@${DOUBLETFINDER_SHA}',
                                 upgrade='never')"

echo "==> Smoke-testing R<->Python anndata bridge"
"$ENV_RSCRIPT" -e 'suppressMessages(library(anndata)); cat("anndata OK\n")'

echo "==> Writing byte-exact lockfile: $HERE/conda.lock.txt"
conda list --name "$ENV_NAME" --explicit --md5 > "$HERE/conda.lock.txt"
echo "    $(wc -l < "$HERE/conda.lock.txt" | tr -d ' ') lines"

# Record the resolved DoubletFinder SHA next to the lockfile so a replay
# can verify both conda and GitHub-sourced components.
echo "DoubletFinder: chris-mcginnis-ucsf/DoubletFinder@$DOUBLETFINDER_SHA" \
    > "$HERE/external_pins.txt"

echo
echo "Setup complete. Activate with:  conda activate $ENV_NAME"
echo "Then run:                       bash run.sh --help"
echo "[setup.sh] Finished $(date '+%Y-%m-%d %H:%M:%S')"
