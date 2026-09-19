#!/usr/bin/env bash
# Compiles every Rust example and runs every Python example as a smoke test.
set -uo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fail=0
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

find "$root/lectures" -mindepth 2 -maxdepth 2 -type d | sort | while read -r dir; do
  name="$(basename "$dir")"
  if [ -f "$dir/example.rs" ]; then
    if ! rustc --edition 2021 -O "$dir/example.rs" -o "$tmp/$name" 2>"$tmp/$name.rs.err"; then
      echo "RUST FAIL: $dir/example.rs"; cat "$tmp/$name.rs.err"; fail=1
    elif ! "$tmp/$name" >/dev/null 2>"$tmp/$name.rs.run.err"; then
      echo "RUST RUNTIME FAIL: $dir/example.rs"; cat "$tmp/$name.rs.run.err"; fail=1
    fi
  fi
  if [ -f "$dir/example.py" ]; then
    if ! python3 "$dir/example.py" >/dev/null 2>"$tmp/$name.py.err"; then
      echo "PYTHON FAIL: $dir/example.py"; cat "$tmp/$name.py.err"; fail=1
    fi
  fi
done

exit $fail
