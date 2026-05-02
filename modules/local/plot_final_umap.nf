process PLOT_FINAL_UMAP {
    tag "final_umap"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'umap_final.png'

    input:
        path h5ad

    output:
        path "umap_final.png", emit: plot

    script:
    """
    python ${projectDir}/bin/plot_final_umap.py \\
        --input ${h5ad} \\
        --plot_output umap_final.png \\
        --color singler_labels
    """
}
