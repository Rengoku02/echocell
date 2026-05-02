process INGEST_H5AD {
    tag "${input.name}"
    publishDir "${params.outdir}/loaded", mode: 'copy', enabled: params.keep_intermediates

    input:
        path input

    output:
        path "complete.h5ad", emit: data

    script:
    def subsample_arg = params.subsample ? "--subsample ${params.subsample}" : ''
    """
    python ${projectDir}/bin/ingest_h5ad_full.py \\
        --h5ad ${input} \\
        --output complete.h5ad \\
        --seed ${params.seed} \\
        ${subsample_arg}
    """
}
