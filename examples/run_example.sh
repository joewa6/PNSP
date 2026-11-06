#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON=${PYTHON:-python}
BUILD_ROOT=${BUILD_ROOT:-"$ROOT/example_output"}

"$PYTHON" "$ROOT/pnsp/build_contaminant.py" \
  "$ROOT/examples/example_contaminants.csv" \
  --output "$BUILD_ROOT/ligands" \
  --microstate-source input

"$PYTHON" "$ROOT/pnsp/build_system.py" \
  --pet "$ROOT/examples/example_pet" \
  --ligand "$BUILD_ROOT/ligands/phenol" \
  --output "$BUILD_ROOT/phenol/rep1" \
  --copies 10 \
  --replica rep1 \
  --seed 314159

"$PYTHON" "$ROOT/pnsp/run_md.py" \
  "$BUILD_ROOT/phenol/rep1" \
  --output "$BUILD_ROOT/phenol/rep1" \
  --production-ns 15

echo "Prepared the example through system construction and AMBER input generation."
echo "Inspect $BUILD_ROOT/phenol/rep1, then run:"
echo "  bash $BUILD_ROOT/phenol/rep1/run_md.sh"
