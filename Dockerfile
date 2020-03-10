FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y gcc g++ make cmake git
RUN apt-get install -y python3-dev python3-venv python3-pip
RUN apt-get install -y libopenmpi-dev libhdf5-openmpi-dev
RUN apt-get install -y freeglut3-dev

ENV HDF5_DIR /usr/lib/x86_64-linux-gnu/hdf5/openmpi
WORKDIR /code/fastpli

CMD git clean -d -f -x && \
   make BUILD=debug install && \
   make test && \
   make examples/requirements && \
   make docs && \
   env/bin/python3 examples/sandbox.py && \
   env/bin/python3 examples/solver.py && \
   env/bin/python3 examples/simpli.py && \
   env/bin/python3 examples/simulation_pipeline.py
