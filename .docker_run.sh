#!/bin/bash

echo "Build"
make BUILD=debug install && \

echo "Test"
make test && \

echo "Docs"
make docs && \

echo "Examples"
make examples/requirements && \
for f in examples/*.py; do echo "$f"; env/bin/python3 $f; done

echo "Readme"
sed -n '/^```sh/,/^```/ p' < README.md | sed '/^```/ d' | source /dev/stdin
