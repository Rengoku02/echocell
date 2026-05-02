process EXPORT_RDS {
    tag "export_rds"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'final.rds'

    input:
        path h5ad

    output:
        path "final.rds", emit: data

    script:
    """
    Rscript ${projectDir}/bin/export_h5ad_to_rds.R \\
        --input ${h5ad} \\
        --output final.rds \\
        --seed ${params.seed}
    """
}
