#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = common \
          C \
          fandc \
          nc_test \
          cdf_format \
          header \
          testcases \
          nonblocking

if BUILD_LARGE_FILE_TEST
SUBDIRS += largefile
endif

if HAS_MPICXX
SUBDIRS += CXX
endif

if HAS_FORTRAN
SUBDIRS += nf_test
if HAVE_MPI_MOD
SUBDIRS += nf90_test
SUBDIRS += F90
endif
endif

if ENABLE_SUBFILING
SUBDIRS += subfile
endif

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
	for i in 2 4 6 8 10 ; do \
	(echo "=============================================")  && \
	(echo "    Parallel testing on $$i MPI processes")  && \
	(echo "=============================================")  && \
	( for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest$$i $$* ; \
	  done ) ; \
	done ;

.PHONY: tests ptest ptests

.NOTPARALLEL:

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

