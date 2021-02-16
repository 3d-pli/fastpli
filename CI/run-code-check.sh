#!/bin/bash
set -uo pipefail

EXIT_STATUS=0

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
   find ./src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' | xargs $CLANGFORMAT --dry-run --Werror
   if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi
else
   echo "WARNING: clang-format not found."
   echo "SKIPPING c/c++ format check."
fi

# Python
env-CI/bin/python3 -m yapf -r -p -d src
if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi
env-CI/bin/python3 -m yapf -r -p -d tests
if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi
env-CI/bin/python3 -m yapf -r -p -d examples
if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi
env-CI/bin/python3 -m flake8
if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi

# Notebooks
find ./examples -iname '*ipynb' | xargs -I {} env-CI/bin/jupyter-nbconvert --clear-output --ClearMetadataPreprocessor.enabled=True --stdout {} | diff {} - &>/dev/null
if [ ! $? -eq 0 ]; then EXIT_STATUS=1; fi

exit $EXIT_STATUS
