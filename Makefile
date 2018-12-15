.DEFAULT: help

.PHONY: help
help:
	@echo "make install"
	@echo "make test"
	@echo "make clean"
	@echo ""
	@echo "make build -- for compile checking"

BUILD := debug

VENV.debug := .venv
VENV.release := .env
VENV := ${VENV.${BUILD}}
PYTHON := ${VENV}/bin/python3
INSTALL.debug := install -v --global-option build --global-option --debug -e .
INSTALL.release := install . -q
INSTALL := ${INSTALL.${BUILD}}

${VENV}:
	python3 -m venv ${VENV}/
	${VENV}/bin/pip3 install --upgrade pip -q
	${VENV}/bin/pip3 install -e . -q

.PHONY: install
install: ${VENV}
	${VENV}/bin/pip3 ${INSTALL}

.PHONY: build
.ONESHELL:
build: clean-build
	mkdir build
	cd build
	cmake ..
	make

.PHONY: test
test:
	${PYTHON} -m unittest discover -s tests -p '*_test.py'

.PHONY: test-docker
test-docker:
	docker build -t fastpli .
	docker run fastpli

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

.PHONY: clean-git
clean-git: check
	git clean src/ -d -f -x -q

.PHONY: check
check:
	@echo -n "Are you sure? [y/N] " && read ans && [ $$ans = y ]
