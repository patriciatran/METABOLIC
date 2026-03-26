process RUN_HMMSEARCH {
    conda 'bioconda::hmmer=3.3.2'
    tag "${meta} - ${hmm_file.baseName}" // Better for troubleshooting

    publishDir "${params.outdir}/intermediate_files/Hmmsearch_Outputs", mode:'copy'

    input:
    // This matches the [meta, faa, hmm] structure from your map
    //tuple val(meta), path(faa), path(hmm_file)
    tuple val(meta), path(faa), path(hmm_file), val(hmm_name)

    output:
    // Emitting a tuple allows you to keep the Sample ID attached to the result
    //tuple val(meta), path("${meta.id}.${hmm_file.baseName}.hmmsearch_result.txt"), emit: results
    path("${meta}.${hmm_name}.hmmsearch_result.txt")

    script:
    // Use .baseName to avoid calling 'basename' in bash
    def hmm_name = hmm_file.baseName
    """
    hmmsearch \
        --tblout ${meta}.${hmm_name}.hmmsearch_result.txt \
        --cpu ${task.cpus} \
        $hmm_file \
        $faa
    """
}