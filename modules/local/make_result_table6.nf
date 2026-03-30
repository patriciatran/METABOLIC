process MAKE_RESULT_TABLE6 {
    label 'process_single'
    publishDir "${params.outdir}/METABOLIC_result_each_spreadsheet", mode: 'copy'

    input:
    path(merops_parsed_results)  // collection of .parsed files from RUN_MEROPS

    output:
    path("METABOLIC_result_worksheet6.tsv"), emit: worksheet6

    script:
    """
    python3 << 'EOF'
import os
import re
from collections import defaultdict

# Parse MEROPS results
merops_hits = defaultdict(lambda: defaultdict(int))      # genome => merops_id => count
merops_hit_names = defaultdict(lambda: defaultdict(str)) # genome => merops_id => hit_names
all_merops_ids = set()
genome_ids = set()

# Process all .parsed files in current directory (collected from RUN_MEROPS)
for filename in sorted(os.listdir('.')):
    if not filename.endswith('.MEROPSout.m8.parsed'):
        continue
    
    # Extract genome ID from filename (e.g., "genome1.MEROPSout.m8.parsed" -> "genome1")
    match = re.match(r'(.+?)\\.MEROPS', filename)
    if not match:
        continue
    
    genome_id = match.group(1)
    genome_ids.add(genome_id)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Parse MEROPS parsed output line (tab-separated: merops_id, count, hit_names)
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            merops_id = parts[0].strip()
            try:
                hit_count = int(parts[1])
            except ValueError:
                continue
            
            hit_names = parts[2].strip() if len(parts) >= 3 else 'None'
            
            all_merops_ids.add(merops_id)
            merops_hits[genome_id][merops_id] = hit_count
            merops_hit_names[genome_id][merops_id] = hit_names

# Sort IDs
sorted_merops_ids = sorted(all_merops_ids)
sorted_genomes = sorted(genome_ids)

# Write worksheet6
with open('METABOLIC_result_worksheet6.tsv', 'w') as out:
    # Header
    header = ['MEROPS peptidase ID']
    for gn_id in sorted_genomes:
        header.append(gn_id + ' Hit numbers')
        header.append(gn_id + ' Hits')
    out.write(chr(9).join(header) + chr(10))
    
    # Body
    for merops_id in sorted_merops_ids:
        row = [merops_id]
        for gn_id in sorted_genomes:
            hit_count = merops_hits[gn_id][merops_id]
            hit_names = merops_hit_names[gn_id][merops_id] if merops_hit_names[gn_id][merops_id] else 'None'
            row.append(str(hit_count))
            row.append(hit_names)
        out.write(chr(9).join(row) + chr(10))

print("Worksheet6 generated with " + str(len(sorted_merops_ids)) + " MEROPS peptidases and " + str(len(sorted_genomes)) + " genomes")
EOF
    """
}
