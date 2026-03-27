process HMM_MOTIF_VALIDATE {
  conda 'python=3.11' 
    tag "${meta} - ${hmm_file.baseName}"

    publishDir "${params.outdir}/intermediate_files/motif_validation", mode: 'copy'

    input:
      tuple val(meta), path(faa), path(hmm_result), path(hmm_file)

    output:
      path("${meta}.${hmm_file.baseName}.motifcheck.tsv")

    script:
    def motif_file = params.motif_file ?: "${params.METABOLIC_db_path}/METABOLIC_template_and_database/motif.txt"
    def motif_pair_file = params.motif_pair_file ?: "${params.METABOLIC_db_path}/METABOLIC_template_and_database/motif.pair.txt"
    """
    set -e
    HMM_NAME=\$(basename ${hmm_file} .hmm)
    OUT="${meta}.\${HMM_NAME}.motifcheck.tsv"

    python - "${faa}" "${hmm_result}" "${motif_file}" "${motif_pair_file}" "\$OUT" <<'PY'
import re
import sys
from pathlib import Path

def load_kv_colon(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            k, v = line.split(':', 1)
            d[k] = v
    return d

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

faa, hmmres, motif_file, motif_pair_file, out_file = sys.argv[1:]
motif = load_kv_colon(motif_file)
motif_pair = load_kv_colon(motif_pair_file)
seqs = read_fasta(faa)

with open(hmmres) as hr, open(out_file, 'w') as out:
    for line in hr:
        if not line or line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 3:
            continue
        hit_id = f[0]
        hmm_name = Path(f[2]).name  # e.g. K00001.hmm
        hmmid = hmm_name[:-4] if hmm_name.endswith('.hmm') else hmm_name
        seq = seqs.get(hit_id)
        if seq is None:
            continue

        status = 'Absent'
        if hmmid in motif:
            pattern = motif[hmmid].replace('X', '[ARNDCQEGHILKMFPSTWYV]')
            status = 'Passed' if re.search(pattern, seq) else 'Failed'
        elif hmmid in motif_pair:
            status = 'Pair-check-needed'

        print(hmmid + "	" + hit_id + "	" + status, file=out)
PY
    """
}