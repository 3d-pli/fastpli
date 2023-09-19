#!/bin/bash
set -euo pipefail

make clean

if [ $# -eq 0 ]; then
   make fastpli
else
   make fastpli $1
fi

python3 -m venv env-CI
env-CI/bin/pip3 install --upgrade pip
env-CI/bin/pip3 install .
env-CI/bin/pip3 install pre-commit
