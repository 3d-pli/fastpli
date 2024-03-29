.PHONY: help
help:
	@echo ----------Installation----------
	@echo make fastpli
	@echo pip3 install .
	@echo -------------Tests--------------
	@echo python3 setup.py test
	@echo -------------Clean--------------
	@echo make clean
	@echo --------------------------------

BUILD := release
VENV := ${if ${venv},${venv},env}
PYTHON3 := python3
PYTHON3_VERSION := $(shell ${PYTHON3} -V | cut -d ' ' -f2)
CMAKE.debug := cmake .. -DCMAKE_BUILD_TYPE=Debug
CMAKE.info := cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
CMAKE.release := cmake .. -DCMAKE_BUILD_TYPE=Release
CMAKE := ${CMAKE.${BUILD}} -DPython3_VERSION=${PYTHON3_VERSION}

MAKE.debug := make
MAKE.info := make -j
MAKE.release := make -j
MAKE := ${MAKE.${BUILD}}

INSTALL.debug := install .
INSTALL.info := install . -q
INSTALL.release := install . -q
INSTALL := ${INSTALL.${BUILD}}

CLANG-FORMAT=clang-format-10

build/:
	mkdir build

build/Makefile: build/
	@if [ ! -f build/Makefile ]; then \
		cd build; \
		echo ${CMAKE}; \
		${CMAKE}; \
	fi

.PHONY: fastpli
fastpli: build/ build/Makefile
	cd build; \
	${MAKE}

.PHONY: clean
clean: clean-build clean-src

################################################################################
############################ DEVELOPER FUNCTIONS ###############################
################################################################################

${VENV}/bin/pip3:
	rm -rf ${VENV}
	${PYTHON3} -m venv ${VENV}
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}/bin/python3:
	rm -rf ${VENV}
	${PYTHON3} -m venv ${VENV}
	${VENV}/bin/pip3 install --upgrade pip -q

.PHONY: ${VENV}
${VENV}: ${VENV}/bin/pip3 ${VENV}/bin/python3

.PHONY: local
local: ${VENV} fastpli
	${VENV}/bin/pip3 ${INSTALL}

.PHONY: development
development: ${VENV} fastpli
	${VENV}/bin/pip3 install -e . -q
	${VENV}/bin/pip3 install yapf pylint flake8 -q
	${VENV}/bin/pip3 install -r examples/requirements.txt -q

.PHONY: uninstall
uninstall:
	@if [ -f ${VENV}/bin/pip3 ]; then \
		if ! ${VENV}/bin/pip3 list | grep -q "fastli"; then \
			${VENV}/bin/pip3 uninstall fastpli -y; \
		fi \
	fi

.PHONY: test
test:
	@# python3 -m unittest tests.analysis.orientation_test.MainTest.test_remap_half_sphere_z
	${VENV}/bin/python3 tests/test.py

.PHONY: h5py-serial
h5py-serial: h5py-clean
	${VENV}/bin/pip3 install h5py

.PHONY: h5py-mpi
h5py-mpi: h5py-clean
	${VENV}/bin/pip3 install -U cython
	HDF5_DIR=${HDF5_DIR} CC=mpicc HDF5_MPI="ON" ${VENV}/bin/pip3 install --no-binary=h5py h5py
	# e.g. HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi

.PHONY: h5py-clean
h5py-clean:
	${VENV}/bin/pip3 uninstall h5py -y

.PHONY: docs
docs: ${VENV} clean-docs
	${VENV}/bin/pip3 -q install -r docs/requirements.txt; \
	cd docs; \
	make html

.PHONY: wiki
wiki: ${VENV}
	${VENV}/bin/pip3 -q install -r wiki/requirements.txt
	rm -rf examples/notebooks
	find examples/ -iname '*.ipynb' | xargs -P 2 -I {} ${VENV}/bin/jupyter-nbconvert --execute --to markdown --output-dir='examples/notebooks' {}
	cd examples/notebooks/; \
	find . -iname '*.md' | xargs -I {} sed -i '2 a ### Note' {} ; \
	find . -iname '*.md' | xargs -I {} sed -i '3G' {} ; \
	find . -iname '*.md' | rev | cut -c4- | rev | cut -c3- | xargs -I {} sed -i '4 a > Jupyter notebook version: [.\/examples\/{}.ipynb](https:\/\/github.com\/3d-pli\/fastpli\/blob\/main\/examples\/{}.ipynb)' {}.md ; \
	find . -iname '*.md' | xargs -I {} sed -i '5G' {} ; \
	find . -iname '*.md' | cut -c3- | xargs -I {} mv {} tutorial-{}

.PHONY: docker
docker:
	./CI/run-docker.sh

.PHONY: format
format:
	find ./src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' | xargs ${CLANG-FORMAT} -i
	${VENV}/bin/pip3 install pre-commit
	pre-commit install
	pre-commit run --all-files

.PHONY: clean-all
clean-all: uninstall clean-build clean-src clean-docs clean-venv

.PHONY: clean-build
clean-build:
	@echo cleaning build
	@rm -rf build
	@rm -f setup.py

.PHONY: clean-docs
clean-docs:
	@echo cleaning docs
	@rm -rf docs/build
	@rm -rf docs/source/_autosummary
	@rm -rf docs/gh-pages

.PHONY: clean-venv
clean-venv:
	@echo cleaning ${VENV}
	@rm -rf ${VENV}
	@rm -rf env-CI

.PHONY: clean-src
clean-src:
	@echo cleaning source
	@rm -f src/fastpli/__version.py
	@rm -f src/include/version.hpp
	@find src/ -name "*.so" -type f | xargs rm -rf
	@find src/ -type f \( -name "*.pyc" -o -name "*.pyo" \) | xargs rm -rf
	@find tests/ -type f \( -name "*.pyc" -o -name "*.pyo" \) | xargs rm -rf
	@find src/ -type d \( -name "__pycache__"  -o -name "*.egg-info" \) | xargs rm -rf
	@find tests/ -type d \( -name "__pycache__" -o -name "*.egg-info" \) | xargs rm -rf
