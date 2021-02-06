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

for f in examples/*.ipynb; do
   echo "converting $f"
   sed '/^\s*"%matplotlib/d' $f | \
   sed  's/plt\.show()/# plt\.show()/g' | \
   env-CI/bin/jupyter-nbconvert --to script --output $f --stdin
   echo "running $f"
   env-CI/bin/python3 $f.py
   rm $f.py
   echo "done $f"
   echo ""
done
