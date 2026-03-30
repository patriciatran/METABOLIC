process MAKE_RESULT_TABLE5 {
    label 'process_single'
    publishDir "${params.outdir}/METABOLIC_result_each_spreadsheet", mode: 'copy'

    input:
    path(dbcan_parsed_results)  // collection of .ps files from RUN_DBCAN

    output:
    path("METABOLIC_result_worksheet5.tsv"), emit: worksheet5

    script:
    """
    python3 << 'EOF'
import os
import re
from collections import defaultdict

# Parse dbCAN results
dbcan_hits = defaultdict(lambda: defaultdict(int))      # genome => cazy_id => count
dbcan_hit_names = defaultdict(lambda: defaultdict(str)) # genome => cazy_id => hit_names
all_cazy_ids = set()
genome_ids = set()

# Process all .ps files in current directory (collected from RUN_DBCAN)
for filename in sorted(os.listdir('.')):
    if not filename.endswith('.dbCAN2.out.dm.ps'):
        continue
    
    # Extract genome ID from filename (e.g., "genome1.dbCAN2.out.dm.ps" -> "genome1")
    match = re.match(r'(.+?)\\.dbCAN', filename)
    if not match:
        continue
    
    genome_id = match.group(1)
    genome_ids.add(genome_id)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse dbCAN parser output line (tab-separated)
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # Check if line starts with CAZy family (GH, PL, CE, AA, CBM, etc.)
            cazy_id_full = parts[0]
            if not any(cazy_id_full.startswith(prefix) for prefix in ['GH', 'PL', 'CE', 'AA', 'CBM']):
                continue
            
            # Extract CAZy ID without .hmm extension
            cazy_match = re.match(r'(\\S+?)\\.hmm', cazy_id_full)
            if cazy_match:
                cazy_id = cazy_match.group(1)
                # Normalize format: ensure 3-digit numbers (e.g., GH1 -> GH001)
                cazy_id_normalized = re.sub(r'([A-Z]+)(\\d+)', lambda m: m.group(1) + m.group(2).zfill(3), cazy_id)
                all_cazy_ids.add(cazy_id_normalized)
                
                # Count hits
                dbcan_hits[genome_id][cazy_id_normalized] += 1
                
                # Store hit names (protein ID is in parts[2])
                if len(parts) >= 3:
                    hit_name = parts[2] if parts[2] else parts[0]
                    if dbcan_hit_names[genome_id][cazy_id_normalized]:
                        dbcan_hit_names[genome_id][cazy_id_normalized] += ';' + hit_name
                    else:
                        dbcan_hit_names[genome_id][cazy_id_normalized] = hit_name

# Sort IDs
sorted_cazy_ids = sorted(all_cazy_ids)
sorted_genomes = sorted(genome_ids)

# Write worksheet5
with open('METABOLIC_result_worksheet5.tsv', 'w') as out:
    # Header
    header = ['CAZyme ID']
    for gn_id in sorted_genomes:
        header.append(gn_id + ' Hit numbers')
        header.append(gn_id + ' Hits')
    out.write(chr(9).join(header) + chr(10))
    
    # Body
    for cazy_id in sorted_cazy_ids:
        row = [cazy_id]
        for gn_id in sorted_genomes:
            hit_count = dbcan_hits[gn_id][cazy_id]
            hit_names = dbcan_hit_names[gn_id][cazy_id] if dbcan_hit_names[gn_id][cazy_id] else 'None'
            row.append(str(hit_count))
            row.append(hit_names)
        out.write(chr(9).join(row) + chr(10))

print("Worksheet5 generated with " + str(len(sorted_cazy_ids)) + " CAZy families and " + str(len(sorted_genomes)) + " genomes")
EOF
    """
}
