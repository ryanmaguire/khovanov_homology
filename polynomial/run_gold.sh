#!/usr/bin/env bash
# One-shot gold regression without make(1).
# Run from repository root:
#   bash polynomial/run_gold.sh
# Or from polynomial/:
#   bash run_gold.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# If script lives in polynomial/, repo root is parent; if copied to root, adjust.
if [[ -d "$SCRIPT_DIR/polynomial" && -f "$SCRIPT_DIR/IntegerMatrix.c" ]]; then
  ROOT="$SCRIPT_DIR"
  POLY="$SCRIPT_DIR/polynomial"
elif [[ -f "$SCRIPT_DIR/../IntegerMatrix.c" ]]; then
  ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
  POLY="$SCRIPT_DIR"
else
  echo "Cannot locate IntegerMatrix.c relative to $SCRIPT_DIR" >&2
  exit 1
fi

CC="${CC:-gcc}"
BUILD="$POLY/build"
mkdir -p "$BUILD"

INC=(-I"$ROOT" -I"$POLY" -I"$POLY/homology" -I"$POLY/polynomial" -I"$POLY/encoder")

SRCS=(
  "$POLY/polynomial/BivariatePoly.c"
  "$POLY/homology/KhPoincare.c"
  "$POLY/homology/TopologyBridge.c"
  "$POLY/encoder/KnotEncoder.c"
  "$ROOT/IntegerMatrix.c"
)

echo "compiling..."
"$CC" -std=c11 -Wall -Wextra -O2 "${INC[@]}" \
  -o "$BUILD/test_gold_khovanov" \
  "$POLY/test/test_gold_khovanov.c" \
  "${SRCS[@]}"

echo "running gold tests..."
"$BUILD/test_gold_khovanov"
echo "OK"
