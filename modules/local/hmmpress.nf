process HMMPRESS_DB {
    // Use bioconda hmmer package (includes hmmpress)
    conda 'bioconda::hmmer=3.3.2'
    tag "${hmm_file.simpleName}"

    publishDir "${output_dir}", mode: 'copy', pattern: '*.h3*'

    input:
    tuple val(output_dir), path(hmm_file)

    output:
    path("${hmm_file}.h3m"), emit: h3m
    path("${hmm_file}.h3i"), emit: h3i
    path("${hmm_file}.h3p"), emit: h3p
    path("${hmm_file}.h3f"), emit: h3f

    script:
    """
    hmmpress ${hmm_file}
    """
}
