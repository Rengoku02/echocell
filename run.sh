#!/usr/bin/env bash
# =============================================================================
#  scRNA-seq analysis pipeline
# =============================================================================
#  Stages: ingest (complete object) -> doublet -> QC -> normalize/HVG/PCA
#          -> pre-integration UMAP -> [pause] -> Harmony? -> UMAP
#          -> SingleR annotation -> final UMAP
#
#  Usage:
#      bash run.sh --from-seurat <rds>  --outdir <out> [options]
#      bash run.sh --from-h5ad  <h5ad>  --outdir <out> [options]
#
#  Flow:
#      --subsample <N>        cell cap during ingest (default: none)
#      --no-prompt            skip the Harmony y/n checkpoint
#      --skip-doublet         bypass DoubletFinder
#      --output-format <fmt>  rds | h5ad | both (default: matches input type)
#      --keep-intermediates   keep per-stage .h5ad files (default: cleanup)
#
#  Doublet:
#      --doublet-rate <F>     expected doublet rate (default: auto from cells)
#
#  QC:
#      --min-genes <N>        min genes per cell
#      --max-genes <N>        max genes per cell
#      --min-umis <N>         min UMIs per cell
#      --max-mito-pct <F>     max mitochondrial %
#      --max-ribo-pct <F>     max ribosomal %
#
#  Embedding / integration:
#      --hvg-flavor <f>       seurat | seurat_v3 | cell_ranger
#      --n-top-genes <N>      number of HVGs
#      --n-pcs <N>            number of principal components
#      --harmony-theta <F>    Harmony diversity penalty
#
#  Annotation:
#      --singler-ref <name>   SingleR reference (e.g. HumanPrimaryCellAtlasData)
#      --singler-labels <c>   reference label column (label.main / label.fine)
#
#  Reproducibility:
#      --seed <N>             random seed for every stochastic step (default 42)
#      --threads <N>          BLAS/OpenMP thread cap (default 1 for determinism)

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
    PYTHON_BIN="$CONDA_PREFIX/bin/python"
else
    PYTHON_BIN="$(command -v python3)"
fi
if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/Rscript" ]]; then
    RSCRIPT_BIN="$CONDA_PREFIX/bin/Rscript"
else
    RSCRIPT_BIN="$(command -v Rscript)"
fi
export RETICULATE_PYTHON="$PYTHON_BIN"

# ----- defaults -----
INPUT_MODE=""
INPUT=""
OUTDIR=""
SUBSAMPLE=""
NO_PROMPT=0
SKIP_DOUBLET=0
HVG_FLAVOR=seurat
OUTPUT_FORMAT=""
KEEP_INTERMEDIATES=0
SEED=42
THREADS=1

# Tunable knobs (empty = use the bin/ script's own default)
DOUBLET_RATE=""
MIN_GENES=""
MAX_GENES=""
MIN_UMIS=""
MAX_MITO_PCT=""
MAX_RIBO_PCT=""
N_TOP_GENES=""
N_PCS=""
HARMONY_THETA=""
SINGLER_REF=""
SINGLER_LABELS=""

need_value() {
    local flag=$1
    local value=${2:-}
    if [[ -z "$value" || "$value" == --* ]]; then
        echo "ERROR: $flag requires a value" >&2
        exit 2
    fi
}

# ----- arg parse -----
while [[ $# -gt 0 ]]; do
    case $1 in
        --from-seurat) need_value "$1" "${2:-}"; INPUT_MODE=seurat; INPUT=$2; shift 2 ;;
        --from-h5ad)   need_value "$1" "${2:-}"; INPUT_MODE=h5ad;   INPUT=$2; shift 2 ;;
        --outdir)      need_value "$1" "${2:-}"; OUTDIR=$2; shift 2 ;;
        --subsample)   need_value "$1" "${2:-}"; SUBSAMPLE=$2; shift 2 ;;
        --no-prompt)   NO_PROMPT=1; shift ;;
        --skip-doublet) SKIP_DOUBLET=1; shift ;;
        --output-format) need_value "$1" "${2:-}"; OUTPUT_FORMAT=$2; shift 2 ;;
        --keep-intermediates) KEEP_INTERMEDIATES=1; shift ;;
        --hvg-flavor)  need_value "$1" "${2:-}"; HVG_FLAVOR=$2; shift 2 ;;
        --doublet-rate)   need_value "$1" "${2:-}"; DOUBLET_RATE=$2; shift 2 ;;
        --min-genes)      need_value "$1" "${2:-}"; MIN_GENES=$2; shift 2 ;;
        --max-genes)      need_value "$1" "${2:-}"; MAX_GENES=$2; shift 2 ;;
        --min-umis)       need_value "$1" "${2:-}"; MIN_UMIS=$2; shift 2 ;;
        --max-mito-pct)   need_value "$1" "${2:-}"; MAX_MITO_PCT=$2; shift 2 ;;
        --max-ribo-pct)   need_value "$1" "${2:-}"; MAX_RIBO_PCT=$2; shift 2 ;;
        --n-top-genes)    need_value "$1" "${2:-}"; N_TOP_GENES=$2; shift 2 ;;
        --n-pcs)          need_value "$1" "${2:-}"; N_PCS=$2; shift 2 ;;
        --harmony-theta)  need_value "$1" "${2:-}"; HARMONY_THETA=$2; shift 2 ;;
        --singler-ref)    need_value "$1" "${2:-}"; SINGLER_REF=$2; shift 2 ;;
        --singler-labels) need_value "$1" "${2:-}"; SINGLER_LABELS=$2; shift 2 ;;
        --seed)           need_value "$1" "${2:-}"; SEED=$2; shift 2 ;;
        --threads)        need_value "$1" "${2:-}"; THREADS=$2; shift 2 ;;
        -h|--help)
            sed -n '3,44p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -z "$INPUT_MODE" || -z "$INPUT" ]] && { echo "ERROR: specify --from-seurat|--from-h5ad"; exit 2; }
[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir required"; exit 2; }

# Default output format matches input type unless explicitly set.
if [[ -z "$OUTPUT_FORMAT" ]]; then
    case $INPUT_MODE in
        seurat) OUTPUT_FORMAT=rds ;;
        h5ad)   OUTPUT_FORMAT=h5ad ;;
    esac
fi
case $OUTPUT_FORMAT in
    rds|h5ad|both) ;;
    *) echo "ERROR: --output-format must be rds|h5ad|both (got '$OUTPUT_FORMAT')" >&2; exit 2 ;;
esac

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Thread caps: same value to every BLAS / OpenMP backend + PythonHashSeed
# keeps dict iteration order stable. Must export BEFORE the stage scripts
# import numpy/scipy.
export OMP_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export VECLIB_MAXIMUM_THREADS=$THREADS
export PYTHONHASHSEED=$SEED

# Mirror all stdout/stderr into log.txt.
LOG="$OUTDIR/log.txt"
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

STARTED_AT="$(date '+%Y-%m-%dT%H:%M:%S%z')"
echo "[run.sh] Started $STARTED_AT"
echo "[run.sh] Mode=$INPUT_MODE  Input=$INPUT  Outdir=$OUTDIR"
echo "[run.sh] Output format: $OUTPUT_FORMAT   Keep intermediates: $KEEP_INTERMEDIATES"
echo "[run.sh] Reproducibility: seed=$SEED threads=$THREADS"
echo "[run.sh] Log: $LOG"

say() { printf "\n\033[1;36m==> %s\033[0m\n" "$*"; }

# ----- Stage 0: ingest complete object -----
mkdir -p "$OUTDIR/loaded" "$OUTDIR/doublet" "$OUTDIR/qc" "$OUTDIR/normalized"
SAMPLE=complete
LOADED="$OUTDIR/loaded/$SAMPLE.h5ad"

case $INPUT_MODE in
    seurat)
        say "Stage 0 — ingest complete Seurat .rds"
        "$RSCRIPT_BIN" "$BIN/ingest_seurat_h5ad.R" --rds "$INPUT" --output "$LOADED" \
            --seed "$SEED" \
            ${SUBSAMPLE:+--subsample $SUBSAMPLE}
        ;;
    h5ad)
        say "Stage 0 — ingest complete .h5ad"
        "$PYTHON_BIN" "$BIN/ingest_h5ad_full.py" --h5ad "$INPUT" --output "$LOADED" \
            --seed "$SEED" \
            ${SUBSAMPLE:+--subsample $SUBSAMPLE}
        ;;
esac

# ----- Stage 1: doublet -> QC -> normalize -----
say "Stage 1 — doublet / QC / normalize"

if [[ $SKIP_DOUBLET -eq 1 ]]; then
    cp "$LOADED" "$OUTDIR/doublet/$SAMPLE.h5ad"
else
    "$RSCRIPT_BIN" "$BIN/run_doubletfinder.R" \
        --input "$LOADED" \
        --output "$OUTDIR/doublet/$SAMPLE.h5ad" \
        --stats_output "$OUTDIR/doublet/${SAMPLE}_stats.txt" \
        --sample_id "$SAMPLE" \
        --seed "$SEED" \
        ${DOUBLET_RATE:+--doublet_rate $DOUBLET_RATE} \
        ${N_PCS:+--pcs $N_PCS}
fi

"$PYTHON_BIN" "$BIN/cell_qc_filter.py" \
    --input "$OUTDIR/doublet/$SAMPLE.h5ad" \
    --output "$OUTDIR/qc/$SAMPLE.h5ad" \
    --stats_output "$OUTDIR/qc/${SAMPLE}_stats.txt" \
    --sample_id "$SAMPLE" \
    ${MIN_GENES:+--min_genes $MIN_GENES} \
    ${MAX_GENES:+--max_genes $MAX_GENES} \
    ${MIN_UMIS:+--min_umis $MIN_UMIS} \
    ${MAX_MITO_PCT:+--max_mito_pct $MAX_MITO_PCT} \
    ${MAX_RIBO_PCT:+--max_ribo_pct $MAX_RIBO_PCT}

"$PYTHON_BIN" "$BIN/normalize_hvg.py" \
    --input "$OUTDIR/qc/$SAMPLE.h5ad" \
    --output "$OUTDIR/normalized/$SAMPLE.h5ad" \
    --sample_id "$SAMPLE" \
    --hvg_flavor "$HVG_FLAVOR" \
    --seed "$SEED" \
    ${N_TOP_GENES:+--n_top_genes $N_TOP_GENES}

# ----- Stage 2: pre-integration UMAP -----
say "Stage 2 — pre-integration UMAP"
"$PYTHON_BIN" "$BIN/umap_checkpoint.py" --inputs "$OUTDIR/normalized/$SAMPLE.h5ad" \
    --output "$OUTDIR/merged_pre.h5ad" \
    --plot_output "$OUTDIR/umap_pre.png" \
    --checkpoint_label pre_integration \
    --seed "$SEED" \
    ${N_TOP_GENES:+--n_top_genes $N_TOP_GENES} \
    ${N_PCS:+--n_pcs $N_PCS}

# ----- Checkpoint: show the UMAP and ask the user -----
if [[ $NO_PROMPT -eq 1 ]]; then
    ans=Y
else
    open "$OUTDIR/umap_pre.png" 2>/dev/null || xdg-open "$OUTDIR/umap_pre.png" 2>/dev/null || true
    printf "\n\033[1;33mPre-integration UMAP: %s\nRun Harmony? [Y/n] \033[0m" "$OUTDIR/umap_pre.png"
    read -r ans
    ans=${ans:-Y}
fi

if [[ "$ans" =~ ^[Yy] ]]; then
    say "Stage 3 — Harmony"
    "$PYTHON_BIN" "$BIN/run_harmony.py" \
        --input "$OUTDIR/merged_pre.h5ad" \
        --output "$OUTDIR/merged_harmony.h5ad" \
        --batch_key sample_id \
        --seed "$SEED" \
        ${N_PCS:+--n_pcs $N_PCS} \
        ${HARMONY_THETA:+--theta $HARMONY_THETA}

    say "Stage 4 — Post-Harmony UMAP"
    "$PYTHON_BIN" "$BIN/post_harmony_umap.py" \
        --input "$OUTDIR/merged_harmony.h5ad" \
        --output "$OUTDIR/merged_umap.h5ad" \
        --plot_output "$OUTDIR/umap_post.png" \
        --seed "$SEED"

    FOR_ANNOTATION="$OUTDIR/merged_umap.h5ad"
else
    say "Skipping Harmony — annotating the uncorrected UMAP"
    FOR_ANNOTATION="$OUTDIR/merged_pre.h5ad"
fi

# ----- Stage 5: SingleR annotation -----
say "Stage 5 — SingleR annotation"
"$RSCRIPT_BIN" "$BIN/run_singler.R" \
    --input "$FOR_ANNOTATION" \
    --output "$OUTDIR/final.h5ad" \
    --seed "$SEED" \
    ${SINGLER_REF:+--ref $SINGLER_REF} \
    ${SINGLER_LABELS:+--labels $SINGLER_LABELS}

say "Stage 6 — final UMAP colored by SingleR labels"
"$PYTHON_BIN" "$BIN/plot_final_umap.py" \
    --input "$OUTDIR/final.h5ad" \
    --plot_output "$OUTDIR/umap_final.png" \
    --color singler_labels

# ----- Stage 7: materialize requested output format(s) -----
say "Stage 7 — emit final output ($OUTPUT_FORMAT)"
case $OUTPUT_FORMAT in
    h5ad)
        # final.h5ad already in place
        ;;
    rds)
        "$RSCRIPT_BIN" "$BIN/export_h5ad_to_rds.R" \
            --input "$OUTDIR/final.h5ad" --output "$OUTDIR/final.rds" --seed "$SEED"
        rm -f "$OUTDIR/final.h5ad"
        ;;
    both)
        "$RSCRIPT_BIN" "$BIN/export_h5ad_to_rds.R" \
            --input "$OUTDIR/final.h5ad" --output "$OUTDIR/final.rds" --seed "$SEED"
        ;;
esac

# ----- Manifest -----
say "Writing run manifest"
PARAMS_JSON="$OUTDIR/.params.json"
cat > "$PARAMS_JSON" <<EOF
{
  "input_mode":   "$INPUT_MODE",
  "input":        "$INPUT",
  "outdir":       "$OUTDIR",
  "subsample":    "${SUBSAMPLE:-}",
  "no_prompt":    $NO_PROMPT,
  "skip_doublet": $SKIP_DOUBLET,
  "output_format": "$OUTPUT_FORMAT",
  "keep_intermediates": $KEEP_INTERMEDIATES,
  "hvg_flavor":   "$HVG_FLAVOR",
  "doublet_rate": "${DOUBLET_RATE:-}",
  "min_genes":    "${MIN_GENES:-}",
  "max_genes":    "${MAX_GENES:-}",
  "min_umis":     "${MIN_UMIS:-}",
  "max_mito_pct": "${MAX_MITO_PCT:-}",
  "max_ribo_pct": "${MAX_RIBO_PCT:-}",
  "n_top_genes":  "${N_TOP_GENES:-}",
  "n_pcs":        "${N_PCS:-}",
  "harmony_theta": "${HARMONY_THETA:-}",
  "singler_ref":  "${SINGLER_REF:-HumanPrimaryCellAtlasData}",
  "singler_labels": "${SINGLER_LABELS:-label.main}",
  "seed":         $SEED,
  "threads":      $THREADS,
  "harmony_ran":  "$ans"
}
EOF

OUTPUT_LIST=""
[[ -f "$OUTDIR/final.h5ad" ]] && OUTPUT_LIST+="$OUTDIR/final.h5ad,"
[[ -f "$OUTDIR/final.rds"  ]] && OUTPUT_LIST+="$OUTDIR/final.rds,"
OUTPUT_LIST+="$OUTDIR/umap_pre.png,$OUTDIR/umap_final.png"
[[ -f "$OUTDIR/umap_post.png" ]] && OUTPUT_LIST+=",$OUTDIR/umap_post.png"

DOUBLET_SHA_FILE="$HERE/external_pins.txt"
DOUBLET_SHA="unknown"
if [[ -f "$DOUBLET_SHA_FILE" ]]; then
    DOUBLET_SHA="$(awk -F@ '/DoubletFinder/ {print $NF}' "$DOUBLET_SHA_FILE")"
fi

"$PYTHON_BIN" "$BIN/write_manifest.py" \
    --outdir "$OUTDIR" \
    --input "$INPUT" \
    --params-json "$PARAMS_JSON" \
    --outputs "$OUTPUT_LIST" \
    --conda-lock "$HERE/conda.lock.txt" \
    --doublet-sha "$DOUBLET_SHA" \
    --started-at "$STARTED_AT"
rm -f "$PARAMS_JSON"

# ----- Cleanup of per-stage intermediates -----
if [[ $KEEP_INTERMEDIATES -eq 0 ]]; then
    say "Cleanup — removing intermediate h5ad files"
    before=$(du -sh "$OUTDIR" 2>/dev/null | awk '{print $1}')
    rm -rf "$OUTDIR/loaded" "$OUTDIR/normalized"
    rm -f "$OUTDIR/doublet/"*.h5ad "$OUTDIR/qc/"*.h5ad
    rm -f "$OUTDIR/merged_pre.h5ad" "$OUTDIR/merged_harmony.h5ad" "$OUTDIR/merged_umap.h5ad"
    after=$(du -sh "$OUTDIR" 2>/dev/null | awk '{print $1}')
    echo "[cleanup] outdir $before -> $after"
else
    echo "[cleanup] --keep-intermediates set; leaving per-stage files in place"
fi

say "Done. Final output(s) in: $OUTDIR"
echo "[run.sh] Finished $(date '+%Y-%m-%d %H:%M:%S')"
