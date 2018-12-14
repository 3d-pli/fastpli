.DEFAULT: help

.PHONY: help
help:
	@echo "make build -- for debugging only"
	@echo "make install"
	@echo "make debug"
	@echo "make clean"

.PHONY: git-modules
git-modules:
ifeq ($(wildcard lib/pybind11/*), )
	git submodule init
	git submodule update
endif

.PHONY: install
install: git-modules
ifeq ($(wildcard .env/.), )
	python3 -m venv .env/
	.env/bin/pip3 install --upgrade pip -q
endif
	.env/bin/pip3 install . -q

.PHONY: debug
debug: git-modules clean-src
ifeq ($(wildcard .venv/.), )
	python3 -m venv .venv/
	.venv/bin/pip3 install --upgrade pip -q
	.venv/bin/pip3 install -r requirements.txt -q
endif
	.venv/bin/pip3 install -v --global-option build --global-option --debug -e .

.PHONY: build
.ONESHELL:
build: git-modules
ifeq ($(wildcard build/.), )
	mkdir build
endif
	cd build
	cmake ..
	make

.PHONY: test
test:
ifneq ($(wildcard .env/.), )
	.env/bin/python3 -m unittest discover -s tests -p '*_test.py'
endif
ifneq ($(wildcard .venv/.), )
	.venv/bin/python3 -m unittest discover -s tests -p '*_test.py'
endif

.PHONY: docker-run
docker-run:
	docker build -t fastpli .
	docker run fastpli

.PHONY: clean-src
clean-src:
	find src/ -name "*egg-info" -exec rm -r {} +
	find src/ -name "*.so" -exec rm {} +
	find src/ -name "__pycache__" -exec rm -r {} +

.PHONY: clean
clean: clean-src
	rm -rf build

.PHONY: clean-all
clean-all: clean-git-check clean
	rm -rf .env
	rm -rf .venv
	git clean -d -f -x -q

.PHONY: clean-git-check
clean-git-check:
	@echo -n "Are you sure? [y/N] " && read ans && [ $$ans = y ]

.PHONY: clean-git
clean-git: clean-git-check
	git clean src/ -d -f -x -q
