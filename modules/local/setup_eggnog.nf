process SETUP_EGGNOG {
    label 'process_single'
    publishDir "${params.METABOLIC_db_path}", mode: 'copy'

    output:
    path("eggnog_data"), emit: eggnog_db

    script:
    """
    set -euo pipefail

    mkdir -p eggnog_data
    cd eggnog_data

    echo "Downloading eggnog-mapper database files..."
    wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
    wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
    wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
    wget http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz
    wget http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz

    echo "Unzip all"
    gunzip eggnog.db.gz
    tar -xf eggnog.taxa.tar.gz
    gunzip eggnog_proteins.dmnd.gz
    tar -xf mmseqs.tar.gz
    tar -xf pfam.tar.gz

    echo "eggnog-mapper database setup complete"
    cd ..
    """
}
