# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

.PHONY: help Makefile
gh-pages:
	rm -rf gh-pages
	make html
	git clone --branch gh-pages --single-branch git@github.com:3d-pli/fastpli.git gh-pages
	rm -r gh-pages/*
	cp -r $(BUILDDIR)/html/* gh-pages/
	touch gh-pages/.nojekyll

# Convert .py to .ipynb (don't seem to need to --execute? perhaps sphinx-build does it...), then build html:
%: Makefile
	@rm -rf $(BUILDDIR)
	@rm -rf _autosummary
	@rm -rf gh-pages
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
