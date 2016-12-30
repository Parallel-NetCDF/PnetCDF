#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = C tutorial
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
TESTS_DIRS = $(SUBDIRS:%=tests-%)
tests: $(TESTS_DIRS)
$(TESTS_DIRS): 
	$(MAKE) $(MFLAGS) -C $(@:tests-%=%) tests

if ENABLE_COVERAGE
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

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

.PHONY: ptest ptests tests

