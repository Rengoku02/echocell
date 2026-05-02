process SINGLER {
    tag "singler"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'final.h5ad'

    input:
        path h5ad

    output:
        path "final.h5ad", emit: data

    script:
    def ref_arg    = params.singler_ref    ? "--ref ${params.singler_ref}"          : ''
    def labels_arg = params.singler_labels ? "--labels ${params.singler_labels}"    : ''
    """
    Rscript ${projectDir}/bin/run_singler.R \\
        --input ${h5ad} \\
        --output final.h5ad \\
        --seed ${params.seed} \\
        ${ref_arg} \\
        ${labels_arg}
    """
}
