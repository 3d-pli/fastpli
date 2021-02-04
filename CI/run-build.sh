#!/bin/bash
set -euo pipefail

make fastpli

env-CI/bin/pip3 install --upgrade pip
env-CI/bin/pip3 install .
