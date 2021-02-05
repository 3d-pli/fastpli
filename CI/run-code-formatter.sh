#!/bin/bash
set -euo pipefail

# C++
# find supported clang-format version
CLANGFORMAT=""
for v in {11..6}; do
   if [ "$v" -lt "7" ]; then
      v="$v.0"
   fi
   if command -v clang-format-$v &>/dev/null; then
      CLANGFORMAT="clang-format-$v"
      break
   fi
done

# run check
if [ ! "$CLANGFORMAT" == "" ]; then
   find ./src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' | xargs $CLANGFORMAT -i
else
   echo "WARNING: clang-format not found."
   echo "SKIPPING c/c++ format check."
fi

# Python
env-CI/bin/python3 -m yapf -i -r -p src
env-CI/bin/python3 -m yapf -i -r -p tests
env-CI/bin/python3 -m yapf -i -r -p examples
env-CI/bin/python3 -m flake8
