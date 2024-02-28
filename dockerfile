FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y gcc g++ cmake make git
RUN apt-get install -y python3-dev python3-venv python3-pip
RUN apt-get install -y libopenmpi-dev freeglut3-dev

RUN pip3 install jupyterlab

RUN git clone https://github.com/3d-pli/fastpli.git /app

WORKDIR /app

RUN make fastpli
RUN pip3 install .
RUN pip3 install -r examples/requirements.txt

CMD ["jupyter", "lab" , "--port=8888", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
