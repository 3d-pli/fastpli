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
make BUILD=release CC=clang-10 CXX=clang++-10 install
make clean
make BUILD=debug CC=gcc-9 CXX=g++-9 install

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
   if [[ $f == *"_mpi.py" ]]; then
      continue
   fi
   if [[ $f == "examples/crossing.py" ]]; then
      continue
   fi
   echo "running $f"
   env/bin/python3 $f
   echo "done $f"
   echo ""
done

# echo "--------------------------------------"
# echo "*************MPI Examples*************"
# echo "--------------------------------------"
# make HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi h5py-mpi
# make examples/requirements
# for f in examples/*_mpi.py; do
#    echo "running $f"
#    mpirun -n 2 env/bin/python3 -m mpi4py $f
#    echo "done $f"
#    echo ""
# done

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
