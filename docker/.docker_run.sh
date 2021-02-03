#!/bin/bash
set -euo pipefail

echo "--------------------------------------"
echo "****************Build*****************"
echo "--------------------------------------"
make fastpli
python3 -m venv env
env/bin/pip3 install --upgrade pip
env/bin/pip3 install .

echo "--------------------------------------"
echo "*****************Test*****************"
echo "--------------------------------------"
make test

echo "--------------------------------------"
echo "***************Examples***************"
echo "--------------------------------------"
env/bin/pip3 install -r examples/requirements.txt
for f in examples/*.py; do
   if [[ $f == *"_mpi.py" ]]; then
      continue
   fi
   if [[ $f == "examples/crossing.py" ]]; then
      continue
   fi
   echo "running $f"
   env/bin/python3 $f
   echo "done $f"
   echo ""
done

echo "--------------------------------------"
echo "****************Format****************"
echo "--------------------------------------"
make clean
env/bin/pip3 install yapf -q
env/bin/pip3 install flake8 -q
make format
if ! git diff --exit-code; then
   exit 1
fi

echo "--------------------------------------"
echo "*****************Done*****************"
echo "--------------------------------------"
