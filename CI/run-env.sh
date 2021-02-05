#!/bin/bash
set -euo pipefail

python3 -m venv env-CI
env-CI/bin/pip3 install yapf -q
env-CI/bin/pip3 install flake8 -q
