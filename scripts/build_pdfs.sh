#!/usr/bin/env bash
# Compiles every lectures/*/*/notes.tex to notes.pdf in place.
set -uo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fail=0

find "$root/lectures" -mindepth 2 -maxdepth 2 -type d | sort | while read -r dir; do
  tex="$dir/notes.tex"
  [ -f "$tex" ] || continue
  echo "==> Building $tex"
  if ! (cd "$dir" && latexmk -xelatex -interaction=nonstopmode -halt-on-error notes.tex >build.log 2>&1); then
    echo "FAILED: $tex (see $dir/build.log)"
    fail=1
  else
    (cd "$dir" && latexmk -c >/dev/null 2>&1)
  fi
done

exit $fail
