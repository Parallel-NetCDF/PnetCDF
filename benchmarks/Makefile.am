#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = C
DIST_SUBDIRS = C FLASH-IO

if HAS_FORTRAN
if HAVE_MPI_MOD
SUBDIRS += FLASH-IO
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

ptest ptests:
	for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest $$* || exit 1 ; \
	done

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	for d in $(DIST_SUBDIRS) ; do \
	   if [ -d $$d/.deps ] ; then rmdir $$d/.deps || true ; fi ; \
	done ; \
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi

.PHONY: ptest ptests tests

