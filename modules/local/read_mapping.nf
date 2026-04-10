process READ_MAPPING {
    tag "${meta.id}"

    conda 'bioconda::bowtie2=2.5.4 bioconda::minimap2=2.28 bioconda::samtools=1.20'
    publishDir "${params.outdir}/intermediate_files/read_mapping", mode: 'copy'

    input:
    tuple val(meta), val(read1), val(read2)
    path gene_files

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: sorted_bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.flagstat.txt"),   emit: flagstat

    script:
    def seq_type = (params.sequencing_type ?: 'illumina').toString()
    def gene_inputs = (gene_files instanceof List) ? gene_files.collect { "\"${it}\"" }.join(' ') : "\"${gene_files}\""

    """
    set -euo pipefail

    cat ${gene_inputs} > All_gene_collections.gene

    if [[ ! -s All_gene_collections.gene ]]; then
        echo "ERROR: All_gene_collections.gene is empty or missing" >&2
        exit 1
    fi

    if [[ "${seq_type}" == "illumina" ]]; then
        bowtie2-build All_gene_collections.gene All_gene_collections.gene.scaffold

        if [[ -n "${read2}" && "${read2}" != "null" ]]; then
            bowtie2 \
                -x All_gene_collections.gene.scaffold \
                -1 "${read1}" \
                -2 "${read2}" \
                -S "${meta.id}.sam" \
                -p ${task.cpus} 
        else
            bowtie2 \
                -x All_gene_collections.gene.scaffold \
                -U "${read1}" \
                -S "${meta.id}.sam" \
                -p ${task.cpus} 
        fi
    else
        preset="map-ont"
        case "${seq_type}" in
            nanopore)      preset="map-ont" ;;
            pacbio)        preset="map-pb" ;;
            pacbio_hifi)   preset="map-hifi" ;;
            pacbio_asm20)  preset="asm20" ;;
            *)             preset="map-ont" ;;
        esac

        minimap2 \
            -ax "\${preset}" \
            -t ${task.cpus} \
            All_gene_collections.gene \
            "${read1}" > "${meta.id}.sam"
    fi

    samtools view -@ ${task.cpus} -bS "${meta.id}.sam" > "${meta.id}.bam"
    samtools sort -@ ${task.cpus} -o "${meta.id}.sorted.bam" "${meta.id}.bam"
    samtools index -@ ${task.cpus} "${meta.id}.sorted.bam"
    samtools flagstat "${meta.id}.sorted.bam" > "${meta.id}.flagstat.txt"
    """
}
