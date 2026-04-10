process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high_memory'

    conda 'bioconda::gtdbtk=2.6.1'
    publishDir "${params.outdir}/intermediate_files/gtdbtk", mode: 'copy'

    input:
    tuple val(meta), val(genome_dir)
    path db
    val gtdbtk_tmp

    output:
    tuple val(meta), path("${prefix}"),                               emit: gtdb_outdir
    tuple val(meta), path("${prefix}/classify/*.summary.tsv"),        emit: summary
    tuple val(meta), path("${prefix}/classify/*.classify.tree"),      emit: tree, optional: true
    tuple val(meta), path("${prefix}/identify/*.markers_summary.tsv"), emit: markers, optional: true
    tuple val(meta), path("${prefix}/align/*.msa.fasta.gz"),          emit: msa, optional: true
    tuple val(meta), path("${prefix}/align/*.user_msa.fasta.gz"),     emit: user_msa, optional: true
    tuple val(meta), path("${prefix}/align/*.filtered.tsv"),          emit: filtered, optional: true
    tuple val(meta), path("${prefix}/identify/*.failed_genomes.tsv"), emit: failed, optional: true
    tuple val(meta), path("${prefix}/${prefix}.log"),                 emit: log
    tuple val(meta), path("${prefix}/${prefix}.warnings.log"),        emit: warnings

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
//def pplacer_scratch = use_pplacer_scratch_dir ? "--scratch_dir pplacer_tmp" : ""
    def tmpdir_arg = gtdbtk_tmp ? "--tmpdir ${gtdbtk_tmp}" : ""   

    """
    set -euo pipefail

    gtdb_data_path="\$(readlink -f "${db}")"

    if [[ ! -d "\${gtdb_data_path}/metadata" ]]; then
        echo "ERROR: --gtdbtk_db_path must point directly to the GTDB release folder containing metadata/." >&2
        echo "Provided path: \${gtdb_data_path}" >&2
        echo "Missing: \${gtdb_data_path}/metadata" >&2
        exit 1
    fi

    export GTDBTK_DATA_PATH="\${gtdb_data_path}"
    echo "Using GTDBTK_DATA_PATH=\${GTDBTK_DATA_PATH}"

    if [[ -n "${gtdbtk_tmp}" ]]; then
        mkdir -p ${gtdbtk_tmp}
    fi

    gtdbtk classify_wf \\
        ${args} \\
        --genome_dir "${genome_dir}" \
        -x fasta \
        --prefix "${prefix}" \\
        --out_dir ${prefix} \\
        --cpus ${task.cpus} \\
        ${tmpdir_arg}

    if [[ -f ${prefix}/gtdbtk.log ]]; then
        mv ${prefix}/gtdbtk.log "${prefix}/${prefix}.log"
    fi

    if [[ -f ${prefix}/gtdbtk.warnings.log ]]; then
        mv ${prefix}/gtdbtk.warnings.log "${prefix}/${prefix}.warnings.log"
    else
        touch "${prefix}/${prefix}.warnings.log"
    fi

    if [[ ! -f "${prefix}/${prefix}.log" ]]; then
        touch "${prefix}/${prefix}.log"
    fi
    """
}
