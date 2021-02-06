# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
DOCS          = gh-pages

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

.PHONY: help Makefile
github:
	rm -rf gh-pages
	git clone --branch gh-pages --single-branch https://github.com/3d-pli/fastpli.git gh-pages
	rm -rf $(BUILDDIR)/*
	make html
	mv $(BUILDDIR)/html/* $(DOCS)/
	touch $(DOCS)/.nojekyll


# Convert .py to .ipynb (don't seem to need to --execute? perhaps sphinx-build does it...), then build html:
%: Makefile
	rm -rf build
	rm -rf source/_autosummary
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
