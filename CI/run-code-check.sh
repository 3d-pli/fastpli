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
      echo "found $CLANGFORMAT"
      break
   fi
done

# run check
if [ ! "$CLANGFORMAT" == "" ]; then
   find ./src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' -exec $CLANGFORMAT -i {} \;
else
   echo "clang-format not found."
   echo "SKIPPING c/c++ format check."
fi

# Python
env-CI/bin/pip3 install yapf -q
env-CI/bin/pip3 install flake8 -q

env-CI/bin/python3 -m yapf -i -r -p --style pep8 src
env-CI/bin/python3 -m yapf -i -r -p --style pep8 tests
env-CI/bin/python3 -m yapf -i -r -p --style pep8 examples
env-CI/bin/python3 -m flake8 --exclude src/fastpli/__version.py src/fastpli
env-CI/bin/python3 -m flake8 examples
env-CI/bin/python3 -m flake8 tests

if ! git diff --exit-code; then
   echo "Check code format failed."
   echo "Please check and commit."
   exit 1
fi
