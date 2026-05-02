process INGEST_SEURAT {
    tag "${input.name}"
    publishDir "${params.outdir}/loaded", mode: 'copy', enabled: params.keep_intermediates

    input:
        path input

    output:
        path "complete.h5ad", emit: data

    script:
    def subsample_arg = params.subsample ? "--subsample ${params.subsample}" : ''
    """
    Rscript ${projectDir}/bin/ingest_seurat_h5ad.R \\
        --rds ${input} \\
        --output complete.h5ad \\
        --seed ${params.seed} \\
        ${subsample_arg}
    """
}
