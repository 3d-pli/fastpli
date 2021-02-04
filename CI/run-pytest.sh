#!/bin/bash
set -euo pipefail

echo "--------------------------------------"
echo "*****************Test*****************"
echo "--------------------------------------"
env-CI/bin/python3 tests/test.py
