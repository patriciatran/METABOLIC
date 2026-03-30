// make_kegg_identifier_results.nf
// Build per-genome KEGG KO summary files:
//   KEGG_identifier_result/<genome>.result.txt
//   KEGG_identifier_result/<genome>.hits.txt

process MAKE_KEGG_IDENTIFIER_RESULTS {
    conda 'python=3.11'

    publishDir "${params.outdir}/KEGG_identifier_result", mode: 'copy'

    input:
    path motifcheck_files
    path faa_files
    path hmm_table_template

    output:
    path "KEGG_identifier_result/*"

    script:
    """
    python - "${hmm_table_template}" "KEGG_identifier_result" <<'PY'
import os
import sys
from collections import defaultdict

hmm_table_template = sys.argv[1]
out_dir = sys.argv[2]

PASSING = {'Absent', 'Passed', 'Pair-check-needed'}
TAB = chr(9)

# ----------------------------
# Inputs from staged files
# ----------------------------
genome_ids = sorted(f[:-4] for f in os.listdir('.') if f.endswith('.faa'))

# hmm_counts[genome][hmm_file] = passing hit count
# hmm_hits[genome][hmm_file]   = comma-joined passing hit ids
hmm_counts = defaultdict(lambda: defaultdict(int))
hmm_hits = defaultdict(lambda: defaultdict(list))
all_hmms = set()

for fname in sorted(os.listdir('.')):
    if not fname.endswith('.motifcheck.tsv'):
        continue

    base = fname[:-len('.motifcheck.tsv')]
    dot = base.rfind('.')
    if dot < 0:
        continue

    genome = base[:dot]
    hmm_name = base[dot + 1:]
    hmm_file = hmm_name + '.hmm'
    all_hmms.add(hmm_file)

    with open(fname) as fh:
        for raw in fh:
            parts = raw.rstrip().split(TAB)
            if len(parts) < 3:
                continue
            hit_id = parts[1]
            status = parts[2]
            if status in PASSING:
                hmm_counts[genome][hmm_file] += 1
                hmm_hits[genome][hmm_file].append(hit_id)


# ----------------------------
# HMM -> KO mapping
# (Perl _get_hmm_2_KO_hash behavior)
# ----------------------------
hmm2ko = {}
with open(hmm_table_template) as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue

        cols = line.split(TAB)
        if len(cols) < 7:
            continue

        hmm_col = cols[5].strip()
        ko_col = cols[6].strip()
        if not hmm_col:
            continue

        if '; ' not in hmm_col:
            if ko_col.startswith('K'):
                hmm2ko[hmm_col] = ko_col + '.hmm'
        else:
            hmms = [x.strip() for x in hmm_col.split('; ') if x.strip()]
            kos = [x.strip() for x in ko_col.split('; ') if x.strip()]
            for h, k in zip(hmms, kos):
                if k.startswith('K'):
                    hmm2ko[h] = k + '.hmm'


# ----------------------------
# Convert HMM-centric summaries to KO-centric summaries
# (matches Perl overwrite behavior when multiple HMMs map to same KO)
# ----------------------------
ko_ids = set()
ko_result = defaultdict(dict)  # genome -> KO -> count (string)
ko_hits = defaultdict(dict)    # genome -> KO -> comma hits

for gn in sorted(genome_ids):
    for hmm in sorted(all_hmms):
        ko_hmm = ''
        if hmm in hmm2ko:
            ko_hmm = hmm2ko[hmm]
        elif hmm.startswith('K') and hmm.endswith('.hmm') and len(hmm) >= 10 and hmm[1:6].isdigit():
            ko_hmm = hmm

        if not ko_hmm:
            continue

        ko_id = ko_hmm[:-4] if ko_hmm.endswith('.hmm') else ko_hmm
        ko_ids.add(ko_id)

        count = hmm_counts[gn].get(hmm, 0)
        hits = ','.join(hmm_hits[gn].get(hmm, []))

        # Perl assigns directly, so later entries overwrite earlier ones
        ko_result[gn][ko_id] = str(count) if count else ''
        ko_hits[gn][ko_id] = hits if hits else ''


# ----------------------------
# Write outputs
# ----------------------------
os.makedirs(out_dir, exist_ok=True)

for gn in sorted(genome_ids):
    out_result = os.path.join(out_dir, f'{gn}.result.txt')
    out_hits = os.path.join(out_dir, f'{gn}.hits.txt')

    with open(out_result, 'w') as f1, open(out_hits, 'w') as f2:
        for ko_id in sorted(ko_ids):
            print(TAB.join([ko_id, ko_result[gn].get(ko_id, '')]), file=f1)
            print(TAB.join([ko_id, ko_hits[gn].get(ko_id, '')]), file=f2)
PY
    """
}
