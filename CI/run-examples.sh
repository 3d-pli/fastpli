#!/bin/bash
set -euo pipefail

env-CI/bin/pip3 install -r examples/requirements.txt
for f in examples/*.py; do
   if [[ $f == *"_mpi.py" ]]; then
      continue
   fi
   if [[ $f == "examples/crossing.py" ]]; then
      continue
   fi
   echo "running $f"
   env-CI/bin/python3 $f
   echo "done $f"
   echo ""
done
