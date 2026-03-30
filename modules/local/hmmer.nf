process RUN_HMMSEARCH {
    conda 'bioconda::hmmer=3.3.2'
    tag "${meta} - ${hmm_name}" // Better for troubleshooting

    publishDir "${params.outdir}/intermediate_files/Hmmsearch_Outputs", mode:'copy'

    input:
    tuple val(meta), path(faa), path(hmm_file), val(hmm_name), val(threshold)

    output:
    tuple val(meta), path(faa), path("${meta}.${hmm_name}.hmmsearch_result.txt"), path(hmm_file), emit: results

    script:
    """
    hmmsearch \
        -T ${threshold} \
        --tblout ${meta}.${hmm_name}.hmmsearch_result.txt \
        --cpu ${task.cpus} \
        $hmm_file \
        $faa
    """
}