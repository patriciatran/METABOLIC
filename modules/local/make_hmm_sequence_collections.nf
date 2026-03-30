// make_hmm_sequence_collections.nf
// Build per-HMM sequence collections (protein + gene) from motifcheck outputs.
// For each HMM, gather all the passing hits across all genomes and write them to a .faa and/or .gene file with headers in the format >genome~~hit_id.
// Output: there should be a collection.faa and a collgction.gene file for each HMM that had at least one passing hit.

process MAKE_HMM_SEQUENCE_COLLECTIONS {
    conda 'python=3.11'

    publishDir "${params.outdir}/Each_HMM_Amino_Acid_Sequence", mode: 'copy'

    input:
    path motifcheck_files
    path faa_files
    path gene_files

    output:
    path "*.collection.faa", optional: true
    path "*.collection.gene", optional: true

    script:
    """
    python - <<'PY'
import os
from collections import defaultdict

PASSING = {'Absent', 'Passed', 'Pair-check-needed'}


def read_fasta(path):
    seqs = {}
    sid = None
    chunks = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sid is not None:
                    seqs[sid] = ''.join(chunks)
                sid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if sid is not None:
        seqs[sid] = ''.join(chunks)
    return seqs


# 1) Build lookup of all sequences with output-style headers: >genome~~hit_id
faa_lookup = {}
gene_lookup = {}

for fname in sorted(os.listdir('.')):
    if fname.endswith('.faa'):
        genome = fname[:-4]
        for hit_id, seq in read_fasta(fname).items():
            faa_lookup[f'>{genome}~~{hit_id}'] = seq
    elif fname.endswith('.gene'):
        genome = fname[:-5]
        for hit_id, seq in read_fasta(fname).items():
            gene_lookup[f'>{genome}~~{hit_id}'] = seq

# 2) Gather accepted hits by HMM from motifcheck files
#    input name format: {genome}.{hmm}.motifcheck.tsv
hmm_hits = defaultdict(set)

for fname in sorted(os.listdir('.')):
    if not fname.endswith('.motifcheck.tsv'):
        continue

    base = fname[:-len('.motifcheck.tsv')]
    dot = base.rfind('.')
    if dot < 0:
        continue

    genome = base[:dot]
    hmm_name = base[dot + 1:]      # e.g. K03686
    hmm_file = hmm_name + '.hmm'   # match Perl naming

    with open(fname) as fh:
        for raw in fh:
            parts = raw.strip().split()
            if len(parts) < 3:
                continue
            hit_id = parts[1]
            status = parts[2]
            if status in PASSING:
                hmm_hits[hmm_file].add(f'>{genome}~~{hit_id}')

# 3) Write one .faa and/or .gene collection file per HMM
for hmm in sorted(hmm_hits):
    headers = sorted(hmm_hits[hmm])

    faa_records = [(h, faa_lookup[h]) for h in headers if h in faa_lookup]
    gene_records = [(h, gene_lookup[h]) for h in headers if h in gene_lookup]

    if faa_records:
        with open(f'{hmm}.collection.faa', 'w') as out:
            for h, seq in faa_records:
                print(h, file=out)
                print(seq, file=out)

    if gene_records:
        with open(f'{hmm}.collection.gene', 'w') as out:
            for h, seq in gene_records:
                print(h, file=out)
                print(seq, file=out)
PY
    """
}
