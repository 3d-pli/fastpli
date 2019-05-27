.PHONY: default
default: install

.PHONY: help
help:
	@echo "make install"
	@echo "make test"
	@echo "make clean"
	@echo ""
	@echo "make build -- for compile checking"

BUILD := release
VENV := env
PYTHON := ${VENV}/bin/python3
INSTALL.debug := install -v --global-option build --global-option --debug -e .
INSTALL.release := install . -q
INSTALL := ${INSTALL.${BUILD}}

${VENV}/bin/pip3: 
	rm -rf ${VENV}
	python3 -m venv ${VENV}/
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}/bin/python3:
	rm -rf ${VENV}
	python3 -m venv ${VENV}/
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}: ${VENV}/bin/pip3 ${VENV}/bin/python3

.PHONY: git-submodules
git-submodules:
	git submodule update --init

.PHONY: requirements
requirements: 
	${VENV}/bin/pip3 install -r requirements.txt -q

.PHONY: install 
install: ${VENV} git-submodules requirements clean-src
	rm -rf build/
	${VENV}/bin/pip3 ${INSTALL}

.PHONY: h5py-serial
h5py-serial:
	${VENV}/bin/pip3 install h5py

.PHONY: h5py-mpi
h5py-mpi:
	${VENV}/bin/pip3 uninstall h5py
	HDF5_DIR=${HDF5_DIR} CC=mpicc HDF5_MPI="ON" ${VENV}/bin/pip3 install --no-binary=h5py h5py
	# e.g. HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi

.PHONY: h5py-clean
h5py-clean:
	${VENV}/bin/pip3 uninstall h5py -y

.ONESHELL:
build/:
	mkdir build
	cd build
	cmake ..

.PHONY: build
.ONESHELL:
build: build/ git-submodules
	cd build
	make

.PHONY: build-j
.ONESHELL:
build-j: build/ git-submodules
	cd build
	make -j

.PHONY: test
test:
	${PYTHON} -m unittest discover -s tests -p '*_test.py'

.PHONY: docker-build
docker-build:
	docker build -t fastpli .
	docker container rm fastpli-test
	docker create --name fastpli-test fastpli

.PHONY: docker
docker: docker-build
	docker cp . fastpli-test:/code/fastpli/
	docker start -i fastpli-test

.PHONY: clean
clean: clean-src clean-build

.PHONY: clean-all
clean-all: clean-git clean-src clean-build clean-venv

.PHONY: clean-build
clean-build:
	rm -rf build

.PHONY: clean-venv
clean-venv:
	rm -rf .env
	rm -rf .venv

.PHONY: clean-src
clean-src:
	find src/ -name "*egg-info" -exec rm -r {} +
	find src/ -name "*.so" -exec rm {} +
	find src/ -name "__pycache__" -exec rm -r {} +
	find tests/ -name "__pycache__" -exec rm -r {} +

.PHONY: clean-git
clean-git: check
	git clean src/ -d -f -x -q

.PHONY: check
check:
	@echo -n "Are you sure? [y/N] " && read ans && [ $$ans = y ]
