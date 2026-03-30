// setup_dbcan.nf
// Prepare dbCAN2 HMM database before RUN_DBCAN.
// If params.dbcan_hmm_db exists, use it; otherwise download from official dbCAN endpoint.

process SETUP_DBCAN {
    conda 'bioconda::hmmer=3.3.2 python=3.11'

    output:
    path 'dbCAN-fam-HMMs.txt*', emit: hmm_db

    script:
    def dbcan_hmm_db = params.dbcan_hmm_db ?: ''
    def dbcan_download_url = 'https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V14.txt'

    """
    set -euo pipefail

    if [[ -n "${dbcan_hmm_db}" && -f "${dbcan_hmm_db}" ]]; then
        cp "${dbcan_hmm_db}" dbCAN-fam-HMMs.txt
    else
        python - "${dbcan_download_url}" <<'PY'
import sys
import urllib.request

url = sys.argv[1]
out = 'dbCAN-fam-HMMs.txt'
urllib.request.urlretrieve(url, out)
print(f'Downloaded {url} -> {out}')

with open(out, 'rb') as fh:
    head = fh.read(128).decode('latin1', 'replace')

if not head.startswith('HMMER'):
    raise SystemExit(
        'Downloaded dbCAN file is not an HMM database (starts with: ' +
        repr(head[:20]) + ')'
    )
PY
    fi

    hmmpress -f dbCAN-fam-HMMs.txt
    """
}
