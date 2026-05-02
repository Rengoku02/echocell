process DOUBLETFINDER {
    tag "complete"
    publishDir "${params.outdir}/doublet", mode: 'copy', pattern: '*_stats.txt'
    publishDir "${params.outdir}/doublet", mode: 'copy', pattern: '*.h5ad', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "complete.h5ad",      emit: data
        path "complete_stats.txt", emit: stats

    script:
    def doublet_rate_arg = params.doublet_rate ? "--doublet_rate ${params.doublet_rate}" : ''
    def n_pcs_arg        = params.n_pcs        ? "--pcs ${params.n_pcs}"               : ''
    """
    Rscript ${projectDir}/bin/run_doubletfinder.R \\
        --input ${h5ad} \\
        --output complete.h5ad \\
        --stats_output complete_stats.txt \\
        --sample_id complete \\
        --seed ${params.seed} \\
        ${doublet_rate_arg} \\
        ${n_pcs_arg}
    """
}
