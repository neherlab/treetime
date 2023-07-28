#!/usr/bin/env bash

set -euo pipefail

cd test

# Remove treetime_examples in case it exists to not fail
rm -rf treetime_examples
git clone https://github.com/neherlab/treetime_examples.git

bash command_line_tests.sh
OUT=$?
if [ "$OUT" != 0 ]; then
  exit 1
fi

pytest test_treetime.py
if [ "$OUT" != 0 ]; then
  exit 1
fi

# Clean up, the 202* is to remove auto-generated output dirs
rm -rf treetime_examples __pycache__ 202*
