#!/bin/bash
set -euo pipefail

make clean
make fastpli $1

env-CI/bin/pip3 install --upgrade pip
env-CI/bin/pip3 install .
