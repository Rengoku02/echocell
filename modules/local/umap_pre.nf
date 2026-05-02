process UMAP_PRE {
    tag "pre_integration"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'umap_pre.png'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merged_pre.h5ad', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "merged_pre.h5ad", emit: data
        path "umap_pre.png",    emit: plot

    script:
    def n_top_genes_arg = params.n_top_genes ? "--n_top_genes ${params.n_top_genes}" : ''
    def n_pcs_arg       = params.n_pcs       ? "--n_pcs ${params.n_pcs}"             : ''
    """
    python ${projectDir}/bin/umap_checkpoint.py \\
        --inputs ${h5ad} \\
        --output merged_pre.h5ad \\
        --plot_output umap_pre.png \\
        --checkpoint_label pre_integration \\
        --seed ${params.seed} \\
        ${n_top_genes_arg} \\
        ${n_pcs_arg}
    """
}
