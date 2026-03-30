process SETUP_MEROPS {
    conda 'bioconda::diamond=2.0 python=3.11'

    label 'process_single'
    publishDir "${params.METABOLIC_db_path}/MEROPS", mode: 'copy'

    output:
    path("pepunit*.dmnd"), emit: merops_db
    path("pepunit.lib"), emit: merops_lib

    script:
    """
    set -euo pipefail

    echo "Downloading MEROPS pepunit.lib database..."
    curl -L --silent \
        "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib" \
        -o pepunit.lib

    # Validate the download
    if [[ ! -s pepunit.lib ]]; then
        echo "ERROR: MEROPS pepunit.lib download failed or file is empty" >&2
        exit 1
    fi

    # Check that it's a valid FASTA file
    if ! head -1 pepunit.lib | grep -q '^>'; then
        echo "ERROR: Downloaded file does not appear to be FASTA format" >&2
        exit 1
    fi

    # Sanitize FASTA to avoid DIAMOND parse errors from unexpected characters
    # (e.g. spaces or other non-amino-acid symbols in sequence lines)
    python3 << 'PY'
import re

inp = 'pepunit.lib'
out = 'pepunit.lib.cleaned'
removed = 0
seq_lines = 0

with open(inp, 'r', encoding='utf-8', errors='replace') as fin, open(out, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            fout.write(line.rstrip(chr(10) + chr(13)) + chr(10))
            continue

        seq = line.strip().upper()
        if not seq:
            continue

        seq_lines += 1
        clean = re.sub(r'[^A-Z]', '', seq)
        removed += (len(seq) - len(clean))

        if clean:
            fout.write(clean + chr(10))

print(f'Sanitized MEROPS FASTA: sequence lines={seq_lines}, removed_characters={removed}')
PY

    mv pepunit.lib.cleaned pepunit.lib

    echo "Creating DIAMOND database from pepunit.lib..."
    diamond makedb \
        --in pepunit.lib \
        -d pepunit \
        -p "${task.cpus}"

    if [[ ! -f pepunit.dmnd ]]; then
        echo "ERROR: DIAMOND database creation failed (pepunit.dmnd not found)" >&2
        ls -la pepunit* >&2
        exit 1
    fi

    echo "MEROPS database setup complete"
    """
}
