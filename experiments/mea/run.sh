#!/bin/sh
# Full MEA-1 perimeter. ~9 min single core, measured. Run fetch.sh first.
set -eu
cd "$(dirname "$0")"; mkdir -p results
for c in 1-1 2-1 3-1 4-1 6-1 8-1; do for d in 10 17 24; do
  python3 mea_analysis.py "data/$c-$d.spk.txt" "$d" "$c" "results/$c-$d.json"
done; done
python3 score.py
