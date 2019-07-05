.PHONY: default
default: install

.PHONY: help
help:
	@echo "make install"
	@echo "make test"
	@echo "make clean"

BUILD := release
VENV := env

CMAKE.debug := cmake .. -DCMAKE_BUILD_TYPE=Debug 
CMAKE.release := cmake .. -DCMAKE_BUILD_TYPE=Release 
CMAKE := ${CMAKE.${BUILD}}

MAKE.debug := make
MAKE.release := make -j
MAKE := ${MAKE.${BUILD}}

INSTALL.debug := install build/.
INSTALL.release := install build/. -q
INSTALL := ${INSTALL.${BUILD}}

${VENV}/bin/pip3: 
	@rm -rf ${VENV}
	python3 -m venv ${VENV}/
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}/bin/python3:
	@rm -rf ${VENV}
	python3 -m venv ${VENV}/
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}: ${VENV}/bin/pip3 ${VENV}/bin/python3

.PHONY: git-submodules
git-submodules:
	git submodule update --init

.PHONY: requirements
requirements: 
	${VENV}/bin/pip3 install -r requirements.txt -q

.PHONY: example/requirements
example/requirements: 
	${VENV}/bin/pip3 install -r example/requirements.txt -q

.PHONY: install 
install: ${VENV} git-submodules build
	${VENV}/bin/pip3 ${INSTALL}

.PHONY: development 
development: ${VENV} git-submodules requirements example/requirements clean-cmake build
	${VENV}/bin/pip3 ${INSTALL}

.ONESHELL:
build/:
	mkdir build

.ONESHELL:
build/Makefile: build/
	@if [ ! -f build/Makefile ]
	then
		cd build
		echo ${CMAKE} 
		${CMAKE}
	fi

.PHONY: build
.ONESHELL:
build: build/ build/Makefile
	cd build
	${MAKE}

.PHONY: test
test:
	${VENV}/bin/python3 -m unittest discover -s tests -p '*_test.py'

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
clean: clean-build clean-venv #clean-src

.PHONY: clean-build
clean-build:
	rm -rf build

.PHONY: clean-venv
clean-venv:
	rm -rf ${VENV}

.PHONY: clean-src
clean-src:
	@echo rm src/**/*.so
	@find src/ -name "*egg-info" -exec rm -r {} +
	@find src/ -name "*.so" -exec rm {} +
	@find src/ -name "__pycache__" -exec rm -r {} +
	@find tests/ -name "__pycache__" -exec rm -r {} +
