process NORMALIZE_HVG {
    tag "complete"
    publishDir "${params.outdir}/normalized", mode: 'copy', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "complete.h5ad", emit: data

    script:
    def n_top_genes_arg = params.n_top_genes ? "--n_top_genes ${params.n_top_genes}" : ''
    """
    python ${projectDir}/bin/normalize_hvg.py \\
        --input ${h5ad} \\
        --output complete.h5ad \\
        --sample_id complete \\
        --hvg_flavor ${params.hvg_flavor} \\
        --seed ${params.seed} \\
        ${n_top_genes_arg}
    """
}
