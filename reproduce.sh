#!/bin/sh
# Rebuild every figure in paper/figures/ from the committed JSON results.
# No simulation is run: the recorded numbers are the artifact.
set -eu
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
if [ ! -d .venv ]; then
  echo "creating .venv"
  "$PY" -m venv .venv
  ./.venv/bin/pip install --quiet --upgrade pip
  ./.venv/bin/pip install --quiet -r requirements.txt
fi
echo "regenerating figures"
./.venv/bin/python paper/make_figures.py
echo
echo "done. To reproduce the tissue analysis as well (~9 min, one core):"
echo "  experiments/mea/fetch.sh   # downloads public recordings, verifies checksums"
echo "  PATH=\$PWD/.venv/bin:\$PATH experiments/mea/run.sh"
