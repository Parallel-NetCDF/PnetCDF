#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

C_SUBDIRS = common \
            C \
            fandc \
            nc_test \
            cdf_format \
            header \
            testcases \
            nonblocking \
            largefile

DIST_SUBDIRS = $(C_SUBDIRS) CXX nf_test nf90_test F90 subfile datawarp

SUBDIRS = $(C_SUBDIRS)

if HAS_MPICXX
SUBDIRS += CXX
endif

if HAS_FORTRAN
SUBDIRS += nf_test
if HAVE_MPI_MOD
SUBDIRS += nf90_test F90
endif
endif

if ENABLE_SUBFILING
SUBDIRS += subfile
endif

if BUILD_DRIVER_DW
SUBDIRS += datawarp
endif

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
	for i in 2 4 6 ; do \
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

.NOTPARALLEL:

