process MANIFEST {
    tag "manifest"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'manifest.json'

    input:
        path input_data
        path output_files
        val started_at
        val doublet_sha

    output:
        path "manifest.json"

    script:
    // Collect every output file path on a single comma-separated string.
    def output_arg = output_files instanceof List ?
        output_files.collect { it.toString() }.join(',') :
        output_files.toString()
    """
    cat > params.json <<'PARAMS'
{
  "input":         "${params.input}",
  "outdir":        "${params.outdir}",
  "subsample":     ${params.subsample == null ? 'null' : "\"${params.subsample}\""},
  "skip_doublet":  ${params.skip_doublet},
  "skip_harmony":  ${params.skip_harmony},
  "output_format": ${params.output_format == null ? 'null' : "\"${params.output_format}\""},
  "keep_intermediates": ${params.keep_intermediates},
  "hvg_flavor":    "${params.hvg_flavor}",
  "doublet_rate":  ${params.doublet_rate == null ? 'null' : "\"${params.doublet_rate}\""},
  "min_genes":     ${params.min_genes == null ? 'null' : "\"${params.min_genes}\""},
  "max_genes":     ${params.max_genes == null ? 'null' : "\"${params.max_genes}\""},
  "min_umis":      ${params.min_umis == null ? 'null' : "\"${params.min_umis}\""},
  "max_mito_pct":  ${params.max_mito_pct == null ? 'null' : "\"${params.max_mito_pct}\""},
  "max_ribo_pct":  ${params.max_ribo_pct == null ? 'null' : "\"${params.max_ribo_pct}\""},
  "n_top_genes":   ${params.n_top_genes == null ? 'null' : "\"${params.n_top_genes}\""},
  "n_pcs":         ${params.n_pcs == null ? 'null' : "\"${params.n_pcs}\""},
  "harmony_theta": ${params.harmony_theta == null ? 'null' : "\"${params.harmony_theta}\""},
  "singler_ref":   ${params.singler_ref == null ? 'null' : "\"${params.singler_ref}\""},
  "singler_labels": ${params.singler_labels == null ? 'null' : "\"${params.singler_labels}\""},
  "seed":          ${params.seed},
  "threads":       ${params.threads}
}
PARAMS

    python ${projectDir}/bin/write_manifest.py \\
        --outdir . \\
        --input ${input_data} \\
        --params-json params.json \\
        --outputs '${output_arg}' \\
        --conda-lock ${projectDir}/conda.lock.txt \\
        --doublet-sha '${doublet_sha}' \\
        --started-at '${started_at}'
    """
}
