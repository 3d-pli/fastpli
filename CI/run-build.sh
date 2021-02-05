#!/bin/bash
set -euo pipefail

make clean

if [ $# -eq 0 ]; then
   make fastpli
else
   make fastpli $1
fi

env-CI/bin/pip3 install --upgrade pip
env-CI/bin/pip3 install .
