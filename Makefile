.DEFAULT: help
.PHONY: help git-modules build install debug clean-src clean clean-all clean-git clean-git-check

help:
	@echo "make build -- for debugging only"
	@echo "make install"
	@echo "make debug"
	@echo "make clean"

git-modules:
ifeq ($(wildcard lib/pybind11/*), )
	git submodule init
	git submodule update
endif

install: git-modules
ifeq ($(wildcard .env/.), )
	python3 -m venv .env/
	.env/bin/pip3 install --upgrade pip -q
endif
	.env/bin/pip3 install . -q

debug: git-modules clean-src
ifeq ($(wildcard .venv/.), )
	python3 -m venv .venv/
	.venv/bin/pip3 install --upgrade pip -q
	.venv/bin/pip3 install -r requirements.txt -q
endif
	.venv/bin/pip3 install -v --global-option build --global-option --debug -e .

.ONESHELL:
build: git-modules
ifeq ($(wildcard build/.), )
	mkdir build
endif
	cd build
	cmake ..
	make

clean-src:
	find src/ -name "*egg-info" -exec rm -r {} +
	find src/ -name "*.so" -exec rm {} +
	find src/ -name "__pycache__" -exec rm -r {} +

clean: clean-src
	rm -rf build

clean-all: clean-git-check clean
	rm -rf .env
	rm -rf .venv
	git clean -d -f -x -q

clean-git-check:
	@echo -n "Are you sure? [y/N] " && read ans && [ $$ans = y ]

clean-git: clean-git-check
	git clean src/ -d -f -x -q