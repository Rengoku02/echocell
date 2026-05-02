process UMAP_POST {
    tag "post_harmony"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'umap_post.png'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merged_umap.h5ad', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "merged_umap.h5ad", emit: data
        path "umap_post.png",    emit: plot

    script:
    """
    python ${projectDir}/bin/post_harmony_umap.py \\
        --input ${h5ad} \\
        --output merged_umap.h5ad \\
        --plot_output umap_post.png \\
        --seed ${params.seed}
    """
}
