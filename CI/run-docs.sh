#!/bin/bash
set -euo pipefail

env-CI/bin/pip3 install -r docs/requirements.txt
source env-CI/bin/activate
cd docs
make html
