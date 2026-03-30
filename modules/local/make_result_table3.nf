// make_result_table3.nf
// Build METABOLIC_result_worksheet3.tsv (KEGG module presence per genome).

process MAKE_RESULT_TABLE3 {
    conda 'python=3.11'

    publishDir "${params.outdir}/METABOLIC_result_each_spreadsheet", mode: 'copy'

    input:
    path motifcheck_files
    path faa_files
    path hmm_table_template
    path ko_module_table
    path ko_module_step_db
    val module_cutoff

    output:
    path "METABOLIC_result_worksheet3.tsv"

    script:
    """
    python - "${hmm_table_template}" "${ko_module_table}" "${ko_module_step_db}" "${module_cutoff}" "METABOLIC_result_worksheet3.tsv" <<'PY'
import os
import re
import sys
from collections import defaultdict

hmm_table_template = sys.argv[1]
ko_module_table = sys.argv[2]
ko_module_step_db = sys.argv[3]
module_cutoff = float(sys.argv[4])
out_file = sys.argv[5]

PASSING = {'Absent', 'Passed', 'Pair-check-needed'}
TAB = chr(9)


# ----------------------------
# Inputs from staged files
# ----------------------------
genome_ids = sorted(f[:-4] for f in os.listdir('.') if f.endswith('.faa'))

# hmm_counts[genome][hmm_file] = hit count
hmm_counts = defaultdict(lambda: defaultdict(int))
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
            parts = raw.strip().split()
            if len(parts) < 3:
                continue
            status = parts[2]
            if status in PASSING:
                hmm_counts[genome][hmm_file] += 1


# ----------------------------
# HMM -> KO mapping
# (Perl _get_hmm_2_KO_hash)
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
# KEGG categories (ko00002.keg)
# ----------------------------
cat2modules = defaultdict(list)
current_cat = ''
with open(ko_module_table) as fh:
    for raw in fh:
        line = raw.rstrip()
        if line.startswith('C    '):
            current_cat = line[5:].strip()
            if current_cat not in cat2modules:
                cat2modules[current_cat] = []
        elif line.startswith('D      '):
            m = re.search(r'(M\\d{5})', line)
            if m and current_cat:
                cat2modules[current_cat].append(m.group(1))


# ----------------------------
# KEGG step DB
# ----------------------------
# kegg_module["M00804+01"] = (k_string, name)
kegg_module = {}
module2step_number = {}
module2name = {}

with open(ko_module_step_db) as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line.startswith('name'):
            continue
        cols = line.split(TAB)
        if len(cols) < 3:
            continue

        name = cols[0]
        k_string = cols[1]
        m_step = cols[2]

        kegg_module[m_step] = (k_string, name)
        m = re.fullmatch(r'(M.+?)\\+(.+?)', m_step)
        if m:
            module = m.group(1)
            step_num = int(m.group(2))
            module2step_number[module] = step_num
            module2name[module] = name


def determine_module_step(k_string, ko_hits):
    ko_set = set(ko_hits)

    # Replace KO ids with 1/0
    expr = re.sub(r'(K\\d+)', lambda m: '1' if m.group(1) in ko_set else '0', k_string)

    # KEGG-style logic approximation:
    # '+' => AND, ',' => OR
    expr = expr.replace('+', ' and ')
    expr = expr.replace(',', ' or ')

    # Keep only safe chars/operators after replacement
    expr = re.sub(r'[^0-1()\\sandor]', ' ', expr)
    expr = re.sub(r'\\s+', ' ', expr).strip()

    if not expr:
        return 0

    try:
        return 1 if eval(expr, {'__builtins__': {}}, {}) else 0
    except Exception:
        return 0


# ----------------------------
# Step presence per genome
# ----------------------------
module_step_result = defaultdict(dict)  # m_step -> genome -> 0/1

for m_step in sorted(kegg_module):
    k_string, _name = kegg_module[m_step]
    for gn in genome_ids:
        ko_hits = []
        for hmm in sorted(all_hmms):
            hmm_new = ''
            if hmm in hmm2ko:
                hmm_new = hmm2ko[hmm]
            elif re.fullmatch(r'K\\d{5}\\.hmm', hmm):
                hmm_new = hmm

            if not hmm_new:
                continue

            hmm_new_wo_ext = hmm_new[:-4] if hmm_new.endswith('.hmm') else hmm_new
            if hmm_counts[gn].get(hmm, 0) > 0:
                ko_hits.append(hmm_new_wo_ext)

        module_step_result[m_step][gn] = determine_module_step(k_string, ko_hits)


# ----------------------------
# Module presence per genome
# ----------------------------
module_result = defaultdict(dict)
for module in sorted(module2step_number):
    total_steps = module2step_number[module]
    for gn in genome_ids:
        present_no = 0
        for m_step in module_step_result:
            if m_step.startswith(module + '+') and module_step_result[m_step].get(gn, 0):
                present_no += module_step_result[m_step][gn]

        ratio = (present_no / total_steps) if total_steps else 0.0
        module_result[module][gn] = 'Present' if ratio >= module_cutoff else 'Absent'


# ----------------------------
# Worksheet 3 output
# ----------------------------
head = ['Module ID', 'Module', 'Module Category']
for gn in genome_ids:
    head.append(f'{gn} Module presence')

with open(out_file, 'w') as out:
    print(TAB.join(head), file=out)

    for module in sorted(module_result):
        cat = ''
        for c in sorted(cat2modules):
            if module in cat2modules[c]:
                cat = c
                break

        row = [module, module2name.get(module, ''), cat]
        for gn in genome_ids:
            row.append(module_result[module].get(gn, 'Absent'))
        print(TAB.join(row), file=out)
PY
    """
}
