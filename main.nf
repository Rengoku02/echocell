#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INGEST_SEURAT   } from './modules/local/ingest_seurat'
include { INGEST_H5AD     } from './modules/local/ingest_h5ad'
include { DOUBLETFINDER   } from './modules/local/doubletfinder'
include { QC_FILTER       } from './modules/local/qc_filter'
include { NORMALIZE_HVG   } from './modules/local/normalize_hvg'
include { UMAP_PRE        } from './modules/local/umap_pre'
include { HARMONY         } from './modules/local/harmony'
include { UMAP_POST       } from './modules/local/umap_post'
include { SINGLER         } from './modules/local/singler'
include { PLOT_FINAL_UMAP } from './modules/local/plot_final_umap'
include { EXPORT_RDS      } from './modules/local/export_rds'
include { MANIFEST        } from './modules/local/manifest'

// Resolve --output-format. Null → match input type.
def resolveOutputFormat() {
    if (params.output_format) return params.output_format
    return params.input.toString().endsWith('.rds') ? 'rds' : 'h5ad'
}

// Read pinned DoubletFinder SHA (best-effort; falls back to 'unknown').
def doubletFinderSha() {
    def f = file("${projectDir}/external_pins.txt")
    if (!f.exists()) return 'unknown'
    def line = f.readLines().find { it.contains('DoubletFinder') }
    if (!line) return 'unknown'
    return line.split('@').last().trim()
}

workflow {
    if (!params.input) {
        error "ERROR: --input is required (Seurat .rds or AnnData .h5ad)"
    }
    def started_at  = java.time.OffsetDateTime.now().toString()
    def out_format  = resolveOutputFormat()
    def doublet_sha = doubletFinderSha()

    log.info "scrna-seq-pipeline"
    log.info "  input         : ${params.input}"
    log.info "  outdir        : ${params.outdir}"
    log.info "  seed          : ${params.seed}"
    log.info "  threads       : ${params.threads}"
    log.info "  skip_doublet  : ${params.skip_doublet}"
    log.info "  skip_harmony  : ${params.skip_harmony}"
    log.info "  output_format : ${out_format}"

    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    // Stage 0 — ingest (branch on file type)
    if (params.input.toString().endsWith('.rds')) {
        ingested = INGEST_SEURAT(input_ch).data
    } else {
        ingested = INGEST_H5AD(input_ch).data
    }

    // Stage 1a — DoubletFinder (or pass-through)
    doublet = params.skip_doublet ? ingested : DOUBLETFINDER(ingested).data

    // Stage 1b — QC
    qc = QC_FILTER(doublet).data

    // Stage 1c — normalize + HVG + PCA
    normalized = NORMALIZE_HVG(qc).data

    // Stage 2 — pre-integration UMAP (always runs)
    umap_pre = UMAP_PRE(normalized)

    // Stage 3+4 — Harmony + post-integration UMAP (or skip)
    if (params.skip_harmony) {
        for_annotation = umap_pre.data
        umap_post_plot = Channel.empty()
    } else {
        harmonized      = HARMONY(umap_pre.data).data
        umap_post       = UMAP_POST(harmonized)
        for_annotation  = umap_post.data
        umap_post_plot  = umap_post.plot
    }

    // Stage 5 — SingleR
    final_h5ad = SINGLER(for_annotation).data

    // Stage 6 — final UMAP plot
    final_plot = PLOT_FINAL_UMAP(final_h5ad).plot

    // Stage 7 — emit requested format(s)
    rds_ch = (out_format in ['rds', 'both']) ? EXPORT_RDS(final_h5ad).data : Channel.empty()

    // Manifest gets the input file plus every output file we wrote
    output_collection = final_plot
        .mix(umap_pre.plot, umap_post_plot, rds_ch)
        .mix(out_format == 'rds' ? Channel.empty() : final_h5ad)
        .collect()

    MANIFEST(input_ch, output_collection, started_at, doublet_sha)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "Done. Outputs in ${params.outdir}"
    } else {
        log.warn "Pipeline finished with errors — see ${params.outdir}/pipeline_trace.txt"
    }
}
