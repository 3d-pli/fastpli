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
VENV := env

CMAKE.debug := cmake .. -DCMAKE_BUILD_TYPE=Debug
CMAKE.info := cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
CMAKE.release := cmake .. -DCMAKE_BUILD_TYPE=Release
CMAKE := ${CMAKE.${BUILD}}

MAKE.debug := make
MAKE.info := make -j
MAKE.release := make -j
MAKE := ${MAKE.${BUILD}}

INSTALL.debug := install .
INSTALL.info := install . -q
INSTALL.release := install . -q
INSTALL := ${INSTALL.${BUILD}}

DOCKER=ubuntu
CLANG-FORMAT=clang-format-9

${VENV}/bin/pip3:
	rm -rf ${VENV}
	python3 -m venv ${VENV}
	${VENV}/bin/pip3 install --upgrade pip -q

${VENV}/bin/python3:
	rm -rf ${VENV}
	python3 -m venv ${VENV}
	${VENV}/bin/pip3 install --upgrade pip -q

.PHONY: ${VENV}
${VENV}: ${VENV}/bin/pip3 ${VENV}/bin/python3

.PHONY: git-submodules
git-submodules:
	git submodule update --init

.PHONY: examples/requirements
examples/requirements:
	${VENV}/bin/pip3 install -r examples/requirements.txt -q

.PHONY: install
install: ${VENV} fastpli
	${VENV}/bin/pip3 ${INSTALL}

.PHONY: development
development: ${VENV} fastpli
	${VENV}/bin/pip3 install -e . -q
	${VENV}/bin/pip3 install yapf -q
	${VENV}/bin/pip3 install pylint -q
	${VENV}/bin/pip3 install -r examples/requirements.txt -q

.PHONY: uninstall
uninstall:
	@if ! ${VENV}/bin/pip3 list | grep -q "fastli"; then \
		${VENV}/bin/pip3 uninstall fastpli -y; \
	fi

build/:
	mkdir build

build/Makefile: build/
	@if [ ! -f build/Makefile ]; then \
		cd build; \
		echo cmake; \
		${CMAKE}; \
	fi

.PHONY: fastpli
fastpli: build/ build/Makefile
	cd build; \
	${MAKE}

.PHONY: test
test:
	${VENV}/bin/python3 setup.py test

.PHONY: h5py-serial
h5py-serial: h5py-clean
	${VENV}/bin/pip3 install h5py

.PHONY: h5py-mpi
h5py-mpi: h5py-clean
	HDF5_DIR=${HDF5_DIR} CC=mpicc HDF5_MPI="ON" ${VENV}/bin/pip3 install --no-binary=h5py h5py
	# e.g. HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi

.PHONY: h5py-clean
h5py-clean:
	${VENV}/bin/pip3 uninstall h5py -y

.PHONY: docker-build
docker-build:
	docker build -t fastpli-${DOCKER} - < docker/${DOCKER}

.PHONY: docker
docker: docker-build
	rm -rf /tmp/fastpli-${DOCKER}
	git clone . /tmp/fastpli-${DOCKER}
	docker stop fastpli-cont-${DOCKER} || true && docker rm fastpli-cont-${DOCKER} || true
	docker create --name fastpli-cont-${DOCKER} fastpli-${DOCKER}
	docker cp /tmp/fastpli-${DOCKER}/. fastpli-cont-${DOCKER}:/code/fastpli
	rm -rf /tmp/fastpli-${DOCKER}
	docker start -i fastpli-cont-${DOCKER}

.PHONY: docker-parallel
docker-parallel: docker-parallel
	@if [ -f /usr/bin/parallel ]; then \
		parallel --halt now,fail=1 'make DOCKER={} docker' ::: archlinux ubuntu; \
	else \
		make DOCKER=archlinux docker; \
		make DOCKER=ubuntu docker; \
	fi

.PHONY: format
format: format-c++ format-py

.PHONY: format-c++
format-c++:
	find src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' -exec ${CLANG-FORMAT} -i {} \; 
	find tests -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' -exec ${CLANG-FORMAT} -i {} \;

.PHONY: format-py
format-py: 
	${VENV}/bin/python3 -m yapf -i -r -p --style google src;
	${VENV}/bin/python3 -m yapf -i -r -p --style google tests;
	${VENV}/bin/python3 -m yapf -i -r -p --style google examples;

# .PHONY: docs
# docs:
# 	${VENV}/bin/pip3 -q install pdoc3
# 	${VENV}/bin/python3 -m pdoc --force --html --output-dir build/docs fastpli
# 	${VENV}/bin/python3 -m pdoc --pdf fastpli | sed '1,/^...$$/d' | \
# 		sed -e 's/{#\([^]]*\)}/<a name="\1"><\/a>/g' | head -n -2 | \
# 		sed 's/  //g' | sed 's/> //g' | sed 's/######/#####/g' > build/docs/fastpli.md

.PHONY: clean
clean: uninstall clean-build clean-src

.PHONY: clean-all
clean-all: clean-build clean-src clean-venv

.PHONY: clean-build
clean-build:
	@echo cleaning build
	@rm -rf build
	@rm -f setup.py

.PHONY: clean-venv
clean-venv:
	@echo cleaning ${VENV}
	@rm -rf ${VENV}

.PHONY: clean-src
clean-src:
	@echo cleaning source
	@rm -f src/version.py
	@rm -f src/include/version.hpp
	@find src/ -name "*.so" -type f | xargs rm -rf
	@find src/ -type f \( -name "*.pyc" -o -name "*.pyo" \) | xargs rm -rf
	@find tests/ -type f \( -name "*.pyc" -o -name "*.pyo" \) | xargs rm -rf
	@find src/ -type d \( -name "__pycache__"  -o -name "*.egg-info" \) | xargs rm -rf
	@find tests/ -type d \( -name "__pycache__" -o -name "*.egg-info" \) | xargs rm -rf
