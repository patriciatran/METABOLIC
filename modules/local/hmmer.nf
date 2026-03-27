process RUN_HMMSEARCH {
    conda 'bioconda::hmmer=3.3.2'
    tag "${meta} - ${hmm_name}" // Better for troubleshooting

    publishDir "${params.outdir}/intermediate_files/Hmmsearch_Outputs", mode:'copy'

    input:
    // This matches the [meta, faa, hmm] structure from your map
    //tuple val(meta), path(faa), path(hmm_file)
    tuple val(meta), path(faa), path(hmm_file), val(hmm_name), val(threshold)

    output:
    // Emitting a tuple allows you to keep the Sample ID attached to the result
    //tuple val(meta), path("${meta.id}.${hmm_file.baseName}.hmmsearch_result.txt"), emit: results
    path("${meta}.${hmm_name}.hmmsearch_result.txt")

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