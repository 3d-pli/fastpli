#!/bin/bash
set -euo pipefail

env-CI/bin/pip3 install pre-commit jupyter nbconvert
env-CI/bin/pip3 install pre-commit jupyter nbconvert
env-CI/bin/pre-commit run --all-files
