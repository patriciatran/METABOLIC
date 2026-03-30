// draw_element_cycles.nf
// Build per-genome R_input pathway files from motifcheck outputs
// and draw nutrient cycling diagrams via draw_biogeochemical_cycles.R.

process DRAW_ELEMENT_CYCLES {
    conda 'python=3.11 conda-forge::r-base=4.3 conda-forge::r-diagram'
    label 'process_single'

    publishDir "${params.outdir}/METABOLIC_Figures", mode: 'copy', saveAs: { it.startsWith('Nutrient_Cycling_Diagrams') ? it : null }
    publishDir "${params.outdir}/METABOLIC_Figures_Input", mode: 'copy', saveAs: { it.startsWith('Nutrient_Cycling_Diagram_Input') ? it : null }

    input:
    path motifcheck_files
    path faa_files
    path r_pathways_file
    path draw_cycles_script

    output:
    path "Nutrient_Cycling_Diagrams", emit: nutrient_cycling_diagrams
    path "Nutrient_Cycling_Diagram_Input", emit: r_input_files

    script:
    """
    set -euo pipefail

    mkdir -p Nutrient_Cycling_Diagram_Input

    python3 - "${r_pathways_file}" <<'PY'
import os
import sys
from collections import defaultdict

r_pathways_file = sys.argv[1]
PASSING = {'Absent', 'Passed', 'Pair-check-needed'}

# Genome IDs from staged FAA files
# *.faa file names are genome IDs in this workflow
genome_ids = sorted(f[:-4] for f in os.listdir('.') if f.endswith('.faa'))

# hmm_counts[genome][hmm_file] = passing hit count
hmm_counts = defaultdict(lambda: defaultdict(int))

for fname in sorted(os.listdir('.')):
    if not fname.endswith('.motifcheck.tsv'):
        continue

    base = fname[:-len('.motifcheck.tsv')]
    dot = base.rfind('.')
    if dot < 0:
        continue

    genome_id = base[:dot]
    hmm_name = base[dot + 1:]
    hmm_file = hmm_name + '.hmm'

    with open(fname) as fh:
        for raw in fh:
            parts = raw.rstrip().split(chr(9))
            if len(parts) >= 3 and parts[2] in PASSING:
                hmm_counts[genome_id][hmm_file] += 1

# Read R pathways table
r_pathways = {}
r_hmm_ids = set()
with open(r_pathways_file) as fh:
    for raw in fh:
        line = raw.strip()
        if not line:
            continue
        cols = line.split(chr(9))
        if len(cols) < 2:
            continue

        step = cols[0]
        hmms = cols[1]
        r_pathways[step] = hmms

        if ';' not in hmms:
            for key in hmms.split(','):
                key = key.strip()
                if key:
                    r_hmm_ids.add(key)
        else:
            for key in hmms.split(';'):
                for key2 in key.split(','):
                    key2 = key2.strip()
                    if key2 and ('NO' not in key2):
                        r_hmm_ids.add(key2)

# Evaluate pathways for each genome
out_dir = 'Nutrient_Cycling_Diagram_Input'
for gn in genome_ids:
    r_input = {}
    for step in sorted(r_pathways.keys()):
        hmms = r_pathways[step]
        present = 0

        if ';' not in hmms:
            for hmm_id in sorted(r_hmm_ids):
                if (hmm_id in hmms) and (hmm_counts[gn].get(hmm_id, 0) > 0):
                    present = 1
                    break
        else:
            hmms_1, hmms_2 = hmms.split(';', 1)

            if 'NO' not in hmms_2:
                logic1 = 0
                logic2 = 0
                for hmm_id in sorted(r_hmm_ids):
                    if (hmm_id in hmms_1) and (hmm_counts[gn].get(hmm_id, 0) > 0):
                        logic1 = 1
                    if (hmm_id in hmms_2) and (hmm_counts[gn].get(hmm_id, 0) > 0):
                        logic2 = 1
                if logic1 and logic2:
                    present = 1
            else:
                logic1 = 0
                logic2 = 1
                for hmm_id in sorted(r_hmm_ids):
                    if (hmm_id in hmms_1) and (hmm_counts[gn].get(hmm_id, 0) > 0):
                        logic1 = 1
                    if (hmm_id in hmms_2) and (hmm_counts[gn].get(hmm_id, 0) > 0):
                        logic2 = 0
                if logic1 and logic2:
                    present = 1

        r_input[step] = present

    out_file = os.path.join(out_dir, f'{gn}.R_input.txt')
    with open(out_file, 'w') as out:
        for step in sorted(r_input.keys()):
            out.write(step + chr(9) + str(r_input[step]) + chr(10))
PY

    Rscript "${draw_cycles_script}" \
            "Nutrient_Cycling_Diagram_Input" \
      "Output" \
      "FALSE" \
      > /dev/null

        mv Output/draw_biogeochem_cycles Nutrient_Cycling_Diagrams
    rm -rf Output
    """
}
