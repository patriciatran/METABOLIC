process RUN_MEROPS {
    conda 'bioconda::diamond=2.0 python=3.11'

    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}/intermediate_files/MEROPS_Files", mode: 'copy'

    input:
    tuple val(meta), path(faa)
    path merops_db_files
    path merops_lib

    output:
    tuple val(meta), path("${meta.id}.MEROPSout.m8"), emit: results
    tuple val(meta), path("${meta.id}.MEROPSout.m8.parsed"), emit: parsed

    script:
    """
    set -euo pipefail

    if [[ ! -f pepunit.dmnd ]]; then
        echo "ERROR: MEROPS DIAMOND database file pepunit.dmnd not found" >&2
        ls -la pepunit* 2>/dev/null || echo "No pepunit files found"
        exit 1
    fi

    if [[ ! -f pepunit.lib ]]; then
        echo "ERROR: MEROPS library file pepunit.lib not found" >&2
        exit 1
    fi

    echo "Running DIAMOND blastp against MEROPS database for ${meta.id}..."
    diamond blastp \
        -d pepunit \
        -q "$faa" \
        -o ${meta.id}.MEROPSout.m8 \
        -k 1 \
        -e 1e-10 \
        --query-cover 80 \
        --id 50 \
        --quiet \
        -p "${task.cpus}"

    echo "Parsing MEROPS results..."
    python3 << 'PYEOF'
import re
import sys

# Build MEROPS ID mapping from library
merops_map = {}
with open('pepunit.lib', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            # Extract MEROPS ID from header (format: >MER_ID description...)
            # The MER_ID is typically enclosed in #...# markers
            mer_match = re.search(r'#([^#]+)#', line)
            if mer_match:
                mer_id = mer_match.group(1)
                merops_map[line[1:].split()[0]] = mer_id

# Parse DIAMOND results and extract MEROPS IDs
parsed_results = {}
with open('${meta.id}.MEROPSout.m8', 'r') as f:
    for line in f:
        parts = line.strip().split(chr(9))
        if len(parts) < 2:
            continue
        
        query_id = parts[0]
        subject_id = parts[1]
        
        # Look up MEROPS ID
        if subject_id in merops_map:
            mer_id = merops_map[subject_id]
            if mer_id not in parsed_results:
                parsed_results[mer_id] = []
            parsed_results[mer_id].append(query_id)

# Write parsed results
with open('${meta.id}.MEROPSout.m8.parsed', 'w') as out:
    for mer_id in sorted(parsed_results.keys()):
        hit_count = len(parsed_results[mer_id])
        hit_names = ';'.join(parsed_results[mer_id])
        out.write(mer_id + chr(9) + str(hit_count) + chr(9) + hit_names + chr(10))

print("Parsed MEROPS results for ${meta.id}", file=sys.stderr)
PYEOF
    """
}
