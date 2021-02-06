#!/bin/bash
set -euo pipefail

echo "--------------------------------------"
echo "****************env-CI****************"
echo "--------------------------------------"
./CI/run-env.sh
echo "-----------------Done-----------------"
echo "--------------------------------------"
echo "****************Build*****************"
echo "--------------------------------------"
if [ $# -eq 0 ]; then
   ./CI/run-build.sh
else
   ./CI/run-build.sh $1
fi
echo "-----------------Done-----------------"
echo "--------------------------------------"
echo "*************Code-Format**************"
echo "--------------------------------------"
./CI/run-code-check.sh
echo "-----------------Done-----------------"
echo "--------------------------------------"
echo "****************Tests*****************"
echo "--------------------------------------"
./CI/run-pytest.sh
echo "-----------------Done-----------------"
echo "--------------------------------------"
echo "***************Examples***************"
echo "--------------------------------------"
./CI/run-examples.sh
echo "-----------------Done-----------------"
