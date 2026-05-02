process HARMONY {
    tag "harmony"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'merged_harmony.h5ad', enabled: params.keep_intermediates

    input:
        path h5ad

    output:
        path "merged_harmony.h5ad", emit: data

    script:
    def n_pcs_arg         = params.n_pcs         ? "--n_pcs ${params.n_pcs}"             : ''
    def harmony_theta_arg = params.harmony_theta ? "--theta ${params.harmony_theta}"     : ''
    """
    python ${projectDir}/bin/run_harmony.py \\
        --input ${h5ad} \\
        --output merged_harmony.h5ad \\
        --batch_key sample_id \\
        --seed ${params.seed} \\
        ${n_pcs_arg} \\
        ${harmony_theta_arg}
    """
}
