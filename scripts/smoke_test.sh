#!/usr/bin/env bash
set -euo pipefail

# Fast end-to-end smoke test (see scripts/smoke_test.py for details).
# Run from repo root:
#   bash scripts/smoke_test.sh
#
# By default, runs --species all via the integrated runner with grid=256
# and validates:
#   - expected outputs exist
#   - validation_summary_machine.csv has required columns
#   - parameters differ across species (guards against defaults regression)

python3 scripts/smoke_test.py --grid 256 --output smoke_results --keep
