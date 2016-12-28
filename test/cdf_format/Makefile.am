#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib -I$(srcdir)/../common
LDADD = $(top_builddir)/src/lib/libpnetcdf.la ../common/libtestutils.la

TESTPROGRAMS = test_inq_format \
               cdf_type \
               dim_cdf12

check_PROGRAMS = $(TESTPROGRAMS)

AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

CLEANFILES = testfile.nc core core.*

EXTRA_DIST = seq_runs.sh test_cdf1.nc test_cdf2.nc test_cdf5.nc

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptests ptest4: $(TESTPROGRAMS)
	$(TESTMPIRUN4) ./test_inq_format $(srcdir)
	$(TESTMPIRUN4) ./cdf_type $(TESTOUTDIR)/testfile.nc
	$(TESTMPIRUN4) ./dim_cdf12 $(TESTOUTDIR)/testfile.nc

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

