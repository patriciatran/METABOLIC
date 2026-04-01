process RUN_EGGNOG {
    conda 'bioconda::eggnog-mapper'

    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}/intermediate_files/eggnog", mode: 'copy'

    input:
    tuple val(meta), path(faa)
    path eggnog_db

    output:
    tuple val(meta), path("${meta.id}.emapper.annotations"),    emit: annotations
    tuple val(meta), path("${meta.id}.emapper.hits"),           emit: hits
    tuple val(meta), path("${meta.id}.emapper.seed_orthologs"), emit: seed_orthologs

    script:
    """
    set -euo pipefail

    echo "Running eggnog-mapper on ${meta.id}..."

    emapper.py \\
        -m diamond \\
        --itype proteins \\
        -i "$faa" \\
        --output "${meta.id}" \\
        --data_dir "$eggnog_db" \\
        --cpu "${task.cpus}" \\
        --override

    echo "eggnog-mapper annotation complete for ${meta.id}"
    """
}
