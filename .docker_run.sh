#!/bin/bash
set -euo pipefail

echo "--------------------------------------"
echo "****************Build*****************"
echo "--------------------------------------"
make BUILD=debug install

# echo "--------------------------------------"
# echo "****************Readme****************"
# echo "--------------------------------------"
# # readme contains installation, examples and mpi
# sed -n '/^```sh/,/^```/ p' < README.md | sed '/^```/ d' | source /dev/stdin
# make h5py-serial # reinstall h5py

echo "Test"
echo "--------------------------------------"
echo "*****************Test*****************"
echo "--------------------------------------"
make test

# echo "--------------------------------------"
# echo "*****************Docs*****************"
# echo "--------------------------------------"
# make docs

echo "--------------------------------------"
echo "***************Examples***************"
echo "--------------------------------------"
make examples/requirements
for f in examples/*.py; do
   echo "running $f"
   echo env/bin/python3 $f
   echo "done $f"
   echo ""
done

echo "--------------------------------------"
echo "*****************Done*****************"
echo "--------------------------------------"
