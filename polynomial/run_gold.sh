#!/usr/bin/env bash
# Run from repo root:  bash polynomial/run_gold.sh
# Or from polynomial/: bash run_gold.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$SCRIPT_DIR/../IntegerMatrix.c" ]]; then
  ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
  POLY="$SCRIPT_DIR"
else
  echo "Cannot locate IntegerMatrix.c next to polynomial/" >&2
  exit 1
fi

CC="${CC:-gcc}"
BUILD="$POLY/build"
mkdir -p "$BUILD"

"$CC" -std=c11 -Wall -Wextra -O2 \
  -I"$ROOT" -I"$POLY" \
  -o "$BUILD/test_gold_khovanov" \
  "$POLY/test_gold_khovanov.c" \
  "$POLY/BivariatePoly.c" \
  "$POLY/KhPoincare.c" \
  "$POLY/TopologyBridge.c" \
  "$POLY/KnotEncoder.c" \
  "$ROOT/IntegerMatrix.c"

echo "running gold tests..."
"$BUILD/test_gold_khovanov"
echo "OK"
