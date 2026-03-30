// make_result_table2.nf
// Build METABOLIC_result_worksheet2.tsv from motifcheck outputs
// using hmm_table_template.txt + hmm_table_template_2.txt logic from METABOLIC-G.pl.
// This makes a table with Presence or Absent labels for each function for all the genomes.
process MAKE_RESULT_TABLE2 {
    conda 'python=3.11'

    publishDir "${params.outdir}/METABOLIC_result_each_spreadsheet", mode: 'copy'

    input:
    path motifcheck_files
    path faa_files
    path hmm_table_template
    path hmm_table_template_2

    output:
    path "METABOLIC_result_worksheet2.tsv"

    script:
    """
    python - "${hmm_table_template}" "${hmm_table_template_2}" "METABOLIC_result_worksheet2.tsv" <<'PY'
import os
import sys
from collections import defaultdict

hmm_table_1 = sys.argv[1]
hmm_table_2 = sys.argv[2]
out_file = sys.argv[3]

PASSING = {'Absent', 'Passed', 'Pair-check-needed'}

# Genome IDs from staged FAA files
# Keep same style as worksheet1 process
genome_ids = sorted(f[:-4] for f in os.listdir('.') if f.endswith('.faa'))

# hmm_counts[genome]["K00001.hmm"] = integer hit count
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
            parts = raw.rstrip().split('\t')
            if len(parts) >= 3 and parts[2] in PASSING:
                hmm_counts[genome_id][hmm_file] += 1

# Parse template1 rows: line_no -> full cols
template1_rows = {}
header1 = []
with open(hmm_table_1) as fh:
    for raw in fh:
        line = raw.rstrip()
        if not line:
            continue
        if line.startswith('#'):
            header1 = line.split('\t')
            continue
        cols = line.split('\t')
        if cols and cols[0].strip():
            template1_rows[cols[0].strip()] = cols

# Parse template2 rows in input order
template2_rows = []
with open(hmm_table_2) as fh:
    for raw in fh:
        line = raw.rstrip()
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        template2_rows.append(cols)

# Header: template1 cols 1..3 (Category, Function, Gene abbreviation)
# then one "Function presence" per genome
head = [header1[i] if i < len(header1) else '' for i in range(1, 4)]
for gn in genome_ids:
    head.append(f"{gn} Function presence")

with open(out_file, 'w') as out:
    print('\t'.join(head), file=out)

    for cols2 in template2_rows:
        # Perl uses cols index [2..4] for first three display columns
        row = [cols2[i] if i < len(cols2) else '' for i in range(2, 5)]

        # line reference(s) into table1 are in column index 1
        # e.g. "123" or "123||124||130"
        line_ref = cols2[1].strip() if len(cols2) > 1 else ''

        for gn in genome_ids:
            presence = 'Absent'

            if not line_ref:
                row.append(presence)
                continue

            if '||' not in line_ref:
                # Single line reference into table1
                cols1 = template1_rows.get(line_ref, [])
                hmm_col = cols1[5].strip() if len(cols1) > 5 else ''

                if hmm_col and ', ' not in hmm_col:
                    if hmm_counts[gn].get(hmm_col, 0) > 0:
                        presence = 'Present'

                elif hmm_col and ', ' in hmm_col:
                    total = 0
                    for hmm in [x.strip() for x in hmm_col.split(', ')]:
                        total += hmm_counts[gn].get(hmm, 0)
                    if total > 0:
                        presence = 'Present'

            else:
                # Multiple line refs: collect all HMMs from those table1 lines,
                # flatten comma-separated HMM cells, then check summed counts.
                hmms = []
                for ref in line_ref.split('||'):
                    ref = ref.strip()
                    if not ref:
                        continue
                    cols1 = template1_rows.get(ref, [])
                    hmm_col = cols1[5].strip() if len(cols1) > 5 else ''

                    if hmm_col and ', ' not in hmm_col:
                        hmms.append(hmm_col)
                    elif hmm_col and ', ' in hmm_col:
                        hmms.extend([x.strip() for x in hmm_col.split(', ')])

                total = 0
                for hmm in hmms:
                    total += hmm_counts[gn].get(hmm, 0)
                if total > 0:
                    presence = 'Present'

            row.append(presence)

        print('\t'.join(row), file=out)
PY
    """
}
