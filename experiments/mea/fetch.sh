#!/bin/sh
# Fetch the MEA-1 perimeter into ./data/ and verify against manifest.tsv.
# The raw recordings are public and are deliberately not vendored in this repo.
set -eu
cd "$(dirname "$0")"; mkdir -p data
grep -v '^#' manifest.tsv | while IFS="$(printf '\t')" read -r c d sha bytes url; do
  [ -n "${c:-}" ] || continue
  f="data/$c-$d.spk.txt"
  if [ ! -f "$f" ]; then
    curl -fsSL "$url" -o "$f.bz2" && bunzip2 -f "$f.bz2"
  fi
  got=$(shasum -a 256 "$f" | cut -d' ' -f1)
  [ "$got" = "$sha" ] || { echo "CHECKSUM MISMATCH $f" >&2; exit 1; }
done
echo "perimeter verified: $(ls data/*.spk.txt | wc -l) recordings"
