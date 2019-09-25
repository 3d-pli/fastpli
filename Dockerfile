FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y gcc g++ make cmake git
RUN apt-get install -y python3-dev python3-venv python3-pip
RUN apt-get install -y libopenmpi-dev freeglut3-dev libhdf5-openmpi-dev

ENV HDF5_DIR /usr/lib/x86_64-linux-gnu/hdf5/openmpi
WORKDIR /code/fastpli

CMD git clean -d -f -x && make BUILD=release install && make tests
