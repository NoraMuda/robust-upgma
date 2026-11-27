#!/usr/bin/env bash
# Reproduce all Python results of the paper.   Usage: bash run_all.sh [workers]
set -euo pipefail
W=${1:-$(python3 -c "import os; print(os.cpu_count())")}
cd "$(dirname "$0")/python/experiments"
python3 simulation.py --reps 100 --workers "$W"
python3 development.py --round 1 --workers "$W"
python3 development.py --round 2 --workers "$W"
python3 benchmarks.py --reps 30 --workers "$W"
python3 null_structure.py
python3 analysis.py
