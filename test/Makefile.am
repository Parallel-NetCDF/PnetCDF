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
            nonblocking

DIST_SUBDIRS = $(C_SUBDIRS) CXX nf_test nf90_test F90 subfile \
               largefile burst_buffer nc4 adios

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

if ENABLE_NETCDF4
   SUBDIRS += nc4
endif

if ENABLE_ADIOS
   SUBDIRS += adios
endif

if ENABLE_BURST_BUFFER
   SUBDIRS += burst_buffer
endif

if RUN_LARGE_FILE_TEST
   SUBDIRS += largefile
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
	for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest $$* || exit 1 ; \
	done

ptests:
	for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptests $$* || exit 1 ; \
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

.PHONY: ptest ptests

.NOTPARALLEL:

