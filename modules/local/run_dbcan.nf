// run_dbcan.nf
// Run dbCAN2 HMM scanning for each genome FAA.

process RUN_DBCAN {
    conda 'bioconda::hmmer=3.3.2 python=3.11'
    tag "$meta.id"

    publishDir "${params.outdir}/intermediate_files/dbCAN2_Files", mode: 'copy'

    input:
    tuple val(meta), path(faa)
    path dbcan_hmm_db_files

    output:
    tuple val(meta), path("${meta.id}.dbCAN2.out.dm"), path("${meta.id}.dbCAN2.out"), path("${meta.id}.dbCAN2.out.dm.ps"), emit: results

    script:
    def dbcan_parser = params.dbcan_parser ?: "${baseDir}/bin/Accessory_scripts/hmmscan-parser-dbCANmeta.py"

    """
    set -euo pipefail

    if [[ ! -f dbCAN-fam-HMMs.txt ]]; then
      echo "ERROR: staged dbCAN database file dbCAN-fam-HMMs.txt not found" >&2
      exit 1
    fi

    if [[ ! -f "${dbcan_parser}" ]]; then
        echo "ERROR: dbCAN parser not found: ${dbcan_parser}" >&2
        echo "Set --dbcan_parser to hmmscan-parser-dbCANmeta.py" >&2
        exit 1
    fi

    hmmscan \
      --domtblout ${meta.id}.dbCAN2.out.dm \
      --cpu ${task.cpus} \
      dbCAN-fam-HMMs.txt \
      "$faa" \
      > ${meta.id}.dbCAN2.out

    python "${dbcan_parser}" ${meta.id}.dbCAN2.out.dm > ${meta.id}.dbCAN2.out.dm.ps
    """
}
