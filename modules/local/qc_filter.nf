process QC_FILTER {
    tag "complete"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*_stats.txt'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.h5ad', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "complete.h5ad",      emit: data
        path "complete_stats.txt", emit: stats

    script:
    def args = []
    if (params.min_genes    != null) args << "--min_genes ${params.min_genes}"
    if (params.max_genes    != null) args << "--max_genes ${params.max_genes}"
    if (params.min_umis     != null) args << "--min_umis ${params.min_umis}"
    if (params.max_mito_pct != null) args << "--max_mito_pct ${params.max_mito_pct}"
    if (params.max_ribo_pct != null) args << "--max_ribo_pct ${params.max_ribo_pct}"
    """
    python ${projectDir}/bin/cell_qc_filter.py \\
        --input ${h5ad} \\
        --output complete.h5ad \\
        --stats_output complete_stats.txt \\
        --sample_id complete \\
        ${args.join(' ')}
    """
}
