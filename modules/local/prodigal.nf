process ANNOTATE_MAGS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h577a1d6_11'
    tag "$meta.id"
    
    // This publishDir handles the "Move to output folder" logic from your Perl script
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(meta), path(fasta)
    path gff2fasta_script

    output:
    tuple val(meta), path("${meta.id}.faa"),  emit: faa
    tuple val(meta), path("${meta.id}.gene"), emit: genes
    tuple val(meta), path("${meta.id}.gff"),  emit: gff

    script:
    """
    # Prodigal uses $task.cpus assigned from the config
    prodigal \\
        -i $fasta \\
        -a ${meta.id}.faa \\
        -o ${meta.id}.gff \\
        -f gff \\
        -p ${params.prodigal_method} \\
        -q

    perl $gff2fasta_script \\
        -g ${meta.id}.gff \\
        -f $fasta \\
        -o ${meta.id}.gene
    """
}
