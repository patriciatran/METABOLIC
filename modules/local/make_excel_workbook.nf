// make_excel_workbook.nf
// Build METABOLIC_result.xlsx from worksheet1-6 TSV files.

process MAKE_EXCEL_WORKBOOK {
    conda 'conda-forge::r-base=4.3 conda-forge::r-openxlsx'
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path worksheet1
    path worksheet2
    path worksheet3
    path worksheet4
    path worksheet5
    path worksheet6
    path excel_script

    output:
    path "METABOLIC_result.xlsx", emit: workbook

    script:
    """
    set -euo pipefail

    Rscript "${excel_script}" "./" "METABOLIC_result.xlsx" > /dev/null
    """
}
