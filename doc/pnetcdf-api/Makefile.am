#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .tex .pdf

AM_DEFAULT_SOURCE_EXT = .tex

API_GUIDE =
if BUILD_DOCS
if HAS_LATEX
    API_GUIDE += pnetcdf-api.pdf
endif
endif

all: $(API_GUIDE)

pnetcdf-api.pdf: pnetcdf-api.tex
	cd $(srcdir) ; $(LATEX) -output-directory=$(abs_builddir) $< && $(LATEX) -output-directory=$(abs_builddir) $<
	$(DVIPDF) pnetcdf-api.dvi $@

EXTRA_DIST = c_api.tex \
             data_mode_api.tex \
             flexible_api.tex \
             latex8.sty \
             Makefile.in \
             multiple_io.tex \
             nonblocking.tex \
             pnetcdf-api.bbl \
             pnetcdf-api.tex

CLEANFILES = *.dvi *.log *.aux *.ps $(API_GUIDE)

# Turn off parallel builds in this directory.
.NOTPARALLEL:

