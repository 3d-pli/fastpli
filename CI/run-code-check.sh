#!/bin/bash
set -euo pipefail

env-CI/bin/pip3 install pre-commit
env-CI/bin/pre-commit run --all-files
