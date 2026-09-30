#!/usr/bin/env zsh
#
# run_batch.sh — Run a list of alphabetical DT codes through the Khovanov
# scanner and write a CSV of results: DT_code, polynomial_string
#
# Usage:
#   ./run_batch.sh                       # use the inline DEFAULT_LIST below
#   ./run_batch.sh -f knots.txt          # read one DT code per line from file
#   ./run_batch.sh -j 4                  # run 4 parallel workers (GNU xargs)
#   ./run_batch.sh -o khovanov.csv       # custom output path
#   ./run_batch.sh -b ./build/FullScanning  # custom binary path
#
# Exit codes per row (encoded in 3rd CSV column if -x set):
#   0  success
#   1  blocked by scanner
#   2  invalid DT code / parse failure
#   3  torsion-regression self-test failure
#   130 SIGINT/ctrl-c
#
# Notes:
#   * Input lines starting with '#' are skipped.
#   * Inline DEFAULT_LIST is read only when -f is NOT set. To paste your own
#     list, either edit DEFAULT_LIST below or use -f somefile.txt.
#   * Default output: ./khovanov_batch.csv
#   * Default binary : ./FullScanning
#   * Optional -j requires GNU xargs (brew install findutils → gxargs);
#     without -j the script runs serially, which is fine for <= ~20 knots.
# ---------------------------------------------------------------------------

set -o pipefail

SCRIPT_DIR="${0:A:h}"
DEFAULT_BIN="${SCRIPT_DIR}/FullScanning"
DEFAULT_OUT="${SCRIPT_DIR}/khovanov_batch.csv"
DEFAULT_ERR_DIR="${SCRIPT_DIR}/_batch_stderr"

BIN="${DEFAULT_BIN}"
OUT="${DEFAULT_OUT}"
ERR_DIR="${DEFAULT_ERR_DIR}"
INPUT_FILE=""
JOBS=1
EXIT_TAG=0

# ---------------------------------------------------------------------------
# Edit DEFAULT_LIST below, or pass -f knots.txt to override.
# ---------------------------------------------------------------------------
DEFAULT_LIST=(
  # one alphabetical DT code per line, e.g.:
  bca
  bfjihgaedc
)

usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -f <file>      Read DT codes from <file> (one per line, '#' comments OK).
                 Without -f the inline DEFAULT_LIST in this script is used.
  -o <file>      Write CSV output to <file>.       Default: ${DEFAULT_OUT:t}
  -b <path>      Path to FullScanning binary.      Default: ./FullScanning
  -e <dir>       Directory for per-entry stderr.   Default: _batch_stderr/
  -j <N>         Parallel workers.                 Default: 1 (serial)
  -x             Append a per-row exit-code column to the CSV.
  -h             Show this help.
EOF
}

while getopts ":f:o:b:e:j:xh" opt; do
  case "${opt}" in
    f) INPUT_FILE="${OPTARG}" ;;
    o) OUT="${OPTARG}" ;;
    b) BIN="${OPTARG}" ;;
    e) ERR_DIR="${OPTARG}" ;;
    j) JOBS="${OPTARG}" ;;
    x) EXIT_TAG=1 ;;
    h) usage ; exit 0 ;;
    :) echo "error: -${OPTARG} requires an argument" >&2 ; usage >&2 ; exit 2 ;;
    \?) echo "error: unknown option -${OPTARG}" >&2 ; usage >&2 ; exit 2 ;;
  esac
done

if [[ ! -x "${BIN}" ]]; then
  echo "error: binary not found or not executable: ${BIN}" >&2
  echo "       run \`clang -std=c11 -O2 ...\` to build FullScanning first." >&2
  exit 2
fi

# Sanity: integer jobs >=1
if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || (( JOBS < 1 )); then
  echo "error: -j N must be a positive integer, got '${JOBS}'" >&2
  exit 2
fi

mkdir -p "${ERR_DIR}"
: > "${OUT}"

# CSV header
if (( EXIT_TAG )); then
  printf 'dt,polynomial,exit_code\n' > "${OUT}"
else
  printf 'dt,polynomial\n' > "${OUT}"
fi

# ---------------------------------------------------------------------------
# Resolve the list of codes into a temp file (one code per line, stripped).
# ---------------------------------------------------------------------------
RAW_LIST="$(mktemp)"
trap 'rm -f "${RAW_LIST}"' EXIT

if [[ -n "${INPUT_FILE}" ]]; then
  if [[ ! -r "${INPUT_FILE}" ]]; then
    echo "error: cannot read input file '${INPUT_FILE}'" >&2
    exit 2
  fi
  # strip comments and blank lines
  sed -E 's/#.*$//g; /^[[:space:]]*$/d' "${INPUT_FILE}" > "${RAW_LIST}"
else
  printf '%s\n' "${DEFAULT_LIST[@]}" \
    | sed -E 's/#.*$//g; /^[[:space:]]*$/d' > "${RAW_LIST}"
fi

TOTAL=$(wc -l < "${RAW_LIST}" | tr -d ' ')
echo "[batch] codes: ${TOTAL}    binary: ${BIN}    output: ${OUT}    jobs: ${JOBS}" >&2

# ---------------------------------------------------------------------------
# Per-code worker. Lives in a helper script so both the serial while-loop
# (in-process) and the xargs parallel workers (separate zsh processes) can
# find it without needing `export -f`, which is flaky across zsh versions.
# ---------------------------------------------------------------------------
WORKER="${SCRIPT_DIR}/.run_batch_worker.zsh"
cat > "${WORKER}" <<'WORKER_EOF'
#!/usr/bin/env zsh
# Arguments: dt_code
# Env:       BIN OUT ERR_DIR EXIT_TAG
set -o pipefail
local dt="$1"
local safe="${dt//\//_}"
local err="${ERR_DIR}/${safe}.err"
local poly
local -i rc

poly="$("${BIN}" --dt "${dt}" 2> "${err}")"
rc=$?

if [[ "${poly}" == *,* ]] || [[ "${poly}" == *\"* ]] || [[ "${poly}" == *$'\n'* ]]; then
  poly="\"${poly//\"/\"\"}\""
fi

if (( EXIT_TAG )); then
  printf '%s,%s,%d\n' "${dt}" "${poly}" "${rc}" >> "${OUT}"
else
  printf '%s,%s\n' "${dt}" "${poly}" >> "${OUT}"
fi

if (( rc != 0 )); then
  echo "  [rc=${rc}] ${dt}   stderr -> ${err}" >&2
fi
exit 0
WORKER_EOF
chmod +x "${WORKER}"
trap 'rm -f "${RAW_LIST}" "${WORKER}"' EXIT

# Serial dispatch inherits these via process environment; parallel dispatch
# passes them explicitly via `env` on the xargs command line.
export BIN OUT ERR_DIR EXIT_TAG

run_one() {
  "${WORKER}" "$1"
}

# ---------------------------------------------------------------------------
# Dispatch. GNU xargs for parallel (-P) if available; otherwise serial.
# ---------------------------------------------------------------------------
if (( JOBS > 1 )); then
  XARGS=""
  if command -v gxargs >/dev/null 2>&1 ; then XARGS=gxargs
  elif xargs --version 2>/dev/null | grep -q "GNU" ; then XARGS=xargs
  fi
  if [[ -z "${XARGS}" ]]; then
    echo "warning: -j ${JOBS} requires GNU xargs (try: brew install findutils)" >&2
    echo "         falling back to serial execution." >&2
    JOBS=1
  fi
fi

if (( JOBS > 1 )); then
  "${XARGS}" -P "${JOBS}" -I {} env \
    BIN="${BIN}" OUT="${OUT}" ERR_DIR="${ERR_DIR}" EXIT_TAG="${EXIT_TAG}" \
    zsh "${WORKER}" {} < "${RAW_LIST}"
else
  while IFS= read -r dt; do
    run_one "${dt}"
  done < "${RAW_LIST}"
fi

echo "[batch] done. see '${OUT}' (${TOTAL} rows)." >&2
