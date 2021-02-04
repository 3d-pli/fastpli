#!/bin/bash
set -euo pipefail

echo "--------------------------------------"
echo "****************Build*****************"
echo "--------------------------------------"
make fastpli

env-CI/bin/pip3 install --upgrade pip
env-CI/bin/pip3 install .
