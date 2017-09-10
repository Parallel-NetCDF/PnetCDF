#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = C tutorial
DIST_SUBDIRS = C tutorial CXX F77 F90

if HAS_MPICXX
SUBDIRS += CXX
endif

if HAS_FORTRAN
SUBDIRS += F77
if HAVE_MPI_MOD
SUBDIRS += F90
endif
endif

EXTRA_DIST = README

# below is the alternative when AM_EXTRA_RECURSIVE_TARGETS is not defined
# TESTS_DIRS = $(SUBDIRS:%=tests-%)
# tests: $(TESTS_DIRS)
# $(TESTS_DIRS): 
# 	$(MAKE) $(MFLAGS) -C $(@:tests-%=%) tests

if BUILD_COVERAGE
   PTEST_SUBDIRS =
else
   PTEST_SUBDIRS = $(SUBDIRS)
endif

ptest:
	(echo "=============================================")  && \
	(echo "    Parallel testing on 4 MPI processes")  && \
	(echo "=============================================")  && \
	( for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest4 $$* ; \
	done ) ;

ptests:
	for i in 3 4 8 ; do \
	(echo "=============================================")  && \
	(echo "    Parallel testing on $$i MPI processes")  && \
	(echo "=============================================")  && \
	( for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest$$i $$* ; \
	done ) ; \
	done ;

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi

.PHONY: ptest ptests

