#!/bin/bash
set -euo pipefail

# echo "--------------------------------------"
# echo "****************Readme****************"
# echo "--------------------------------------"
# if cat /etc/os-release | grep -iFq ubuntu; then
#    DELETE="pacman"
# elif cat /etc/os-release | grep -iFq archlinux; then
#    DELETE="apt"
# fi

# sed -n '/^```sh/,/^```/ p' < README.md | sed '/^```/ d' | sed "/$DELETE/d" | source /dev/stdin
# make clean

echo "--------------------------------------"
echo "****************Build*****************"
echo "--------------------------------------"
make BUILD=release CC=clang-9 CXX=clang++-9 install
make clean
make BUILD=debug CC=gcc-8 CXX=g++-8 install

echo "--------------------------------------"
echo "*****************Test*****************"
echo "--------------------------------------"
make test

# echo "--------------------------------------"
# echo "*****************Docs*****************"
# echo "--------------------------------------"
# make docs

echo "--------------------------------------"
echo "***************Examples***************"
echo "--------------------------------------"
make examples/requirements
for f in examples/*.py; do
   echo "running $f"
   echo env/bin/python3 $f
   echo "done $f"
   echo ""
done

echo "--------------------------------------"
echo "****************Format****************"
echo "--------------------------------------"
make clean
env/bin/pip3 install yapf -q
make format
if ! git diff --exit-code; then
   exit 1
fi

echo "--------------------------------------"
echo "*****************Done*****************"
echo "--------------------------------------"
