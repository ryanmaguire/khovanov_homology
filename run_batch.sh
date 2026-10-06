#!/usr/bin/env zsh
# Run FullScanning once per alphabetical DT code in a text file.
# Blank lines, whitespace, and # comments are ignored; duplicates are retained.
# Results: dt,polynomial,exit_code (the polynomial is the integral free part).
# Failed knots are recorded and do not prevent the remaining knots from running.
#
# Example:
#   zsh run_batch.sh -f knots.txt -b ./FullScanning.current -o results.csv -j 4
#
# Batch exit status: 0 all succeeded; 1 scanner failures; 2 setup/I/O failures.
# Interrupts exit with 130 (INT) or 143 (TERM). No Makefile is required.

emulate -R zsh
setopt PIPE_FAIL

SCRIPT_DIR="${0:A:h}"
BIN="${SCRIPT_DIR}/FullScanning"
OUT="${SCRIPT_DIR}/khovanov_batch.csv"
ERR_DIR="${SCRIPT_DIR}/_batch_stderr"
INPUT_FILE=""
JOBS=1
FORCE=0
WORK=""
STAGING=""
DISPATCH_PID=""

usage() {
  cat <<EOF
Usage: zsh $0 -f knots.txt [OPTIONS]

  -f file   Required: one alphabetical DT code per line; # comments allowed.
  -b path   Scanner executable. Default: ${BIN}
  -o file   Output CSV. Default: ${OUT}
  -e dir    Parent directory for per-run stderr logs. Default: ${ERR_DIR}
  -j N      Parallel workers. Default: 1.
  -F        Explicitly replace an existing output CSV.
  -x        Accepted for compatibility; exit_code is always included.
  -h        Show help.
EOF
}

die() { print -ru2 -- "error: $*"; exit 2; }

cleanup() {
  # Each active worker traps TERM and stops its scanner child.
  if [[ -n "$WORK" && -d "$WORK" ]]; then
    local pid_file worker_pid
    for pid_file in "$WORK"/*.pid(N); do
      worker_pid=$(<"$pid_file")
      [[ "$worker_pid" == <-> ]] && kill -TERM "$worker_pid" 2>/dev/null
    done
  fi
  if [[ -n "$DISPATCH_PID" ]]; then
    kill -TERM "$DISPATCH_PID" 2>/dev/null
    wait "$DISPATCH_PID" 2>/dev/null
  fi
  [[ -n "$STAGING" ]] && rm -f -- "$STAGING"
  [[ -n "$WORK" ]] && rm -rf -- "$WORK"
  return 0
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

while getopts ':f:o:b:e:j:xFh' opt; do
  case "$opt" in
    f) INPUT_FILE="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    b) BIN="$OPTARG" ;;
    e) ERR_DIR="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    F) FORCE=1 ;;
    x) ;;
    h) usage; exit 0 ;;
    :) die "-$OPTARG requires a value" ;;
    \?) die "unknown option -$OPTARG; use -h" ;;
  esac
done
shift $((OPTIND - 1))
(( $# == 0 )) || die "unexpected arguments; supply input with -f"
[[ -n "$INPUT_FILE" ]] || die "an input file is required: -f knots.txt"
[[ -f "$INPUT_FILE" && -r "$INPUT_FILE" ]] || die "cannot read input file: $INPUT_FILE"
[[ -f "$BIN" && -x "$BIN" ]] || die "scanner is not executable: $BIN"
[[ "$JOBS" =~ ^[0-9]+$ ]] && (( ${#JOBS} <= 10 && JOBS >= 1 && JOBS <= 2147483647 )) \
  || die "-j must be a positive integer no greater than 2147483647"

INPUT_FILE="${INPUT_FILE:A}"
BIN="${BIN:A}"
OUT="${OUT:A}"
ERR_DIR="${ERR_DIR:A}"
[[ "$OUT" != "$INPUT_FILE" && "$OUT" != "$BIN" ]] \
  || die "output must not replace the input file or scanner executable"
[[ ! -d "$OUT" ]] || die "output is a directory: $OUT"
if (( ! FORCE )) && [[ -e "$OUT" || -L "$OUT" ]]; then
  die "output already exists: $OUT (choose a new filename or use -F)"
fi

# Validate and normalize all input before creating or replacing output.
WORK=$(mktemp -d "${TMPDIR:-/tmp}/kh-batch.XXXXXX") || die "cannot create temporary directory"
TASKS="$WORK/tasks.txt"
: > "$TASKS" || die "cannot create task list"
typeset -i LINE_NUMBER=0 TOTAL=0
while IFS= read -r line || [[ -n "$line" ]]; do
  (( ++LINE_NUMBER ))
  line="${line%%\#*}"
  dt="${line//[[:space:]]/}"
  [[ -n "$dt" ]] || continue
  [[ "$dt" =~ ^[A-Za-z]+$ ]] || die "line $LINE_NUMBER is not an alphabetical DT code"
  (( ++TOTAL ))
  printf '%d %s\n' "$TOTAL" "$dt" >> "$TASKS" || die "cannot write task list"
done < "$INPUT_FILE"
(( TOTAL > 0 )) || die "input contains no DT codes"

mkdir -p -- "${OUT:h}" "$ERR_DIR" || die "cannot create output/log directories"
RUN_ERR=$(mktemp -d "$ERR_DIR/run.XXXXXX") || die "cannot create stderr directory"
STAGING=$(mktemp "${OUT:h}/.${OUT:t}.XXXXXX") || die "cannot create output file"
printf 'dt,polynomial,exit_code\n' > "$STAGING" || die "cannot write CSV header"

# Workers use separate files. Only the parent writes the final CSV.
WORKER="$WORK/worker.zsh"
cat > "$WORKER" <<'WORKER_EOF'
#!/usr/bin/env zsh
emulate -R zsh
bin="$1"; work="$2"; errors="$3"; id="$4"; dt="$5"
scanner_pid=""
cleanup_worker() {
  if [[ -n "$scanner_pid" ]]; then
    kill -TERM "$scanner_pid" 2>/dev/null
    wait "$scanner_pid" 2>/dev/null
  fi
  rm -f -- "$work/$id.pid"
}
trap cleanup_worker EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP
print -r -- "$$" > "$work/$id.pid" || exit 2
printf -v log_name '%06d_%s.err' "$id" "$dt"
err="$errors/$log_name"
"$bin" --coeff Z --format polynomial --dt "$dt" \
  > "$work/$id.stdout" 2> "$err" &
scanner_pid=$!
wait "$scanner_pid"
rc=$?
scanner_pid=""
poly=""
if (( rc == 0 )); then
  poly=$(<"$work/$id.stdout")
  if [[ -z "$poly" || "$poly" == *$'\n'* || "$poly" == *$'\r'* ]]; then
    print -r -- 'Batch error: scanner returned no single-line polynomial.' >> "$err" || exit 2
    poly=""
    rc=125
  fi
fi
# Failed computations deliberately have an empty polynomial field.
if [[ "$poly" == *,* || "$poly" == *\"* ]]; then
  poly="\"${poly//\"/\"\"}\""
fi
printf '%s,%s,%d\n' "$dt" "$poly" "$rc" > "$work/$id.csv" || exit 2
printf '%d\n' "$rc" > "$work/$id.rc" || exit 2
exit 0
WORKER_EOF
(( $? == 0 )) || die "cannot write worker script"

print -ru2 -- "[batch] codes: $TOTAL; jobs: $JOBS; binary: $BIN"
print -ru2 -- "[batch] stderr logs: $RUN_ERR"
typeset -i DISPATCH_RC=0
if (( JOBS > 1 )); then
  # Both macOS xargs and GNU xargs support -P and -n.
  command -v xargs >/dev/null || die "parallel execution requires xargs"
  xargs -P "$JOBS" -n 2 zsh "$WORKER" "$BIN" "$WORK" "$RUN_ERR" < "$TASKS" &
  DISPATCH_PID=$!
  wait "$DISPATCH_PID"
  DISPATCH_RC=$?
  DISPATCH_PID=""
else
  while read -r i dt; do
    zsh "$WORKER" "$BIN" "$WORK" "$RUN_ERR" "$i" "$dt" &
    DISPATCH_PID=$!
    wait "$DISPATCH_PID" || DISPATCH_RC=1
    DISPATCH_PID=""
  done < "$TASKS"
fi

typeset -i SUCCEEDED=0 FAILED=0 INFRA_FAILURE=0
(( DISPATCH_RC != 0 )) && INFRA_FAILURE=1
while read -r i dt; do
  rc=""
  [[ -r "$WORK/$i.rc" ]] && rc=$(<"$WORK/$i.rc")
  if [[ "$rc" != <-> || ! -r "$WORK/$i.csv" ]]; then
    printf '%s,,125\n' "$dt" >> "$STAGING" || die "cannot write CSV"
    print -ru2 -- "[batch] missing worker result: entry $i ($dt)"
    INFRA_FAILURE=1
    (( ++FAILED ))
    continue
  fi
  cat "$WORK/$i.csv" >> "$STAGING" || die "cannot write CSV"
  if (( rc == 0 )); then
    (( ++SUCCEEDED ))
  else
    (( ++FAILED ))
    (( rc == 125 )) && INFRA_FAILURE=1
    print -ru2 -- "[batch] failed: entry $i ($dt), exit $rc"
  fi
done < "$TASKS"

# Staging is on the same filesystem as OUT. A hard link provides atomic
# no-overwrite publication, including when two batches choose the same output.
if (( FORCE )); then
  mv -f -- "$STAGING" "$OUT" || die "cannot publish CSV: $OUT"
else
  ln -- "$STAGING" "$OUT" || die "cannot publish CSV; output may now exist: $OUT"
  rm -f -- "$STAGING" || die "cannot remove staging file"
fi
STAGING=""
print -ru2 -- "[batch] finished: $SUCCEEDED succeeded, $FAILED failed; CSV: $OUT"
(( INFRA_FAILURE )) && exit 2
(( FAILED )) && exit 1
exit 0
