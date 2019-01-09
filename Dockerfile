FROM ubuntu

RUN apt-get update
RUN apt-get install -y gcc g++ make cmake git
RUN apt-get install -y python3-dev python3-venv python3-pip
RUN apt-get install -y freeglut3-dev

WORKDIR /code/fastpli

CMD git clean -d -f -x && make git-submodules && make BUILD=release install && make BUILD=release test
