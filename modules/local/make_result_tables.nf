// make_result_tables.nf
// Aggregates per-genome, per-HMM motifcheck results and the HMM template
// to produce METABOLIC_result_worksheet1.tsv (mirrors METABOLIC-G.pl logic).

process MAKE_RESULT_TABLES {
    conda 'python=3.11'

    publishDir "${params.outdir}/METABOLIC_result_each_spreadsheet", mode: 'copy'

    input:
    path motifcheck_files   // collected list – every *.motifcheck.tsv
    path faa_files          // collected list – every *.faa  (used to enumerate genomes)
    path hmm_table_template // hmm_table_template.txt

    output:
    path "METABOLIC_result_worksheet1.tsv"

    script:
    """
    python - "${hmm_table_template}" "METABOLIC_result_worksheet1.tsv" <<'PY'
import sys
import os
from collections import defaultdict

template_file = sys.argv[1]
out_file      = sys.argv[2]

# Statuses from HMM_MOTIF_VALIDATE that represent an accepted hit
# 'Absent'           -> no motif check was needed; hit always counts
# 'Passed'           -> motif regex matched; hit counts
# 'Pair-check-needed'-> pending pair-HMM check; counted conservatively
PASSING = {'Absent', 'Passed', 'Pair-check-needed'}

# ── 1. Enumerate all genome IDs from FAA files staged into work dir ───────────
genome_ids = sorted(f[:-4] for f in os.listdir('.') if f.endswith('.faa'))

# ── 2. Aggregate motifcheck results ───────────────────────────────────────────
# File naming: {genome_id}.{HMM_NAME}.motifcheck.tsv
# Columns:     hmmid   hit_protein_id   status
# hmm_hits[genome_id][hmm_file_name] -> list of passing protein IDs
hmm_hits = defaultdict(lambda: defaultdict(list))

for fname in sorted(os.listdir('.')):
    if not fname.endswith('.motifcheck.tsv'):
        continue
    base = fname[:-len('.motifcheck.tsv')]
    dot  = base.rfind('.')
    if dot < 0:
        continue
    genome_id = base[:dot]
    hmm_name  = base[dot + 1:]
    hmm_file  = hmm_name + '.hmm'          # e.g. K00001.hmm

    with open(fname) as fh:
        for raw in fh:
            parts = raw.strip().split('\t')
            if len(parts) >= 3 and parts[2] in PASSING:
                hmm_hits[genome_id][hmm_file].append(parts[1])

# ── 3. Parse hmm_table_template.txt ───────────────────────────────────────────
# Column layout (0-based after tab-split):
#   0  #Entry   1  Category   2  Function   3  Gene abbreviation
#   4  Gene name   5  Hmm file   6  Corresponding KO   7  Reaction
#   8  Substrate   9  Product   10  Hmm detecting threshold
header_cols   = []
template_rows = {}

with open(template_file) as fh:
    for raw in fh:
        line = raw.strip()
        if not line:
            continue
        if line.startswith('#'):
            header_cols = line.split('\t')
        else:
            cols  = line.split('\t')
            entry = cols[0].strip() if cols else ''
            if entry:
                template_rows[entry] = cols

# ── 4. Write worksheet 1 ──────────────────────────────────────────────────────
# Header: template columns 1-10 (Category -> Hmm detecting threshold),
#         then three columns per genome: Hmm presence / Hit numbers / Hits
head = [header_cols[i] if i < len(header_cols) else '' for i in range(1, 11)]
for gn in genome_ids:
    head.append(gn + ' Hmm presence')
    head.append(gn + ' Hit numbers')
    head.append(gn + ' Hits')

with open(out_file, 'w') as out:
    print('\t'.join(head), file=out)

    for entry in sorted(template_rows.keys()):
        cols    = template_rows[entry]
        row     = [cols[i] if i < len(cols) else '' for i in range(1, 11)]
        hmm_col = cols[5].strip() if len(cols) > 5 else ''

        for gn in genome_ids:
            if not hmm_col:
                # Template row has no associated HMM
                row.extend(['Absent', '0', ''])
                continue

            if ', ' not in hmm_col:
                # ── Single HMM entry (e.g. "K00001.hmm") ──
                h = hmm_hits[gn].get(hmm_col, [])
                row.append('Present' if h else 'Absent')
                row.append(str(len(h)))
                row.append(';'.join(h) if h else 'None')

            else:
                # ── Multiple HMMs, comma-space separated ──
                # e.g. "K00823.hmm, K07250.hmm, K13524.hmm"
                all_hits  = []
                real_hits = 0
                for hmm in [x.strip() for x in hmm_col.split(', ')]:
                    h = hmm_hits[gn].get(hmm, [])
                    if h:
                        real_hits += len(h)
                        all_hits.extend(h)
                    else:
                        all_hits.append('None')
                row.append('Present' if real_hits else 'Absent')
                row.append(str(real_hits))
                row.append(';'.join(all_hits))

        print('\t'.join(row), file=out)
PY
    """
}
