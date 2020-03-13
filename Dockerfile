FROM archlinux/base

RUN pacman --noconfirm --needed -Syu gcc make cmake git
RUN pacman --noconfirm --needed -Syu python python-pipenv python-pip
RUN pacman --noconfirm --needed -Syu openmpi hdf5-openmpi
RUN pacman --noconfirm --needed -Syu freeglut
# glu

WORKDIR /code/fastpli

CMD [ "./.docker_run.sh" ]