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

TESTPROGRAMS = header_consistency

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

# consistency check should be run on more than one process
TESTS = seq_runs.sh

EXTRA_DIST = seq_runs.sh

CLEANFILES = *.nc core core.*

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptest4: $(PROGS)
	set -e ; for i in $(PROGS); do ( \
	$(TEST_MPIRUN_4) ./$$i $(TEST_OUTDIR)/testfile.nc ; \
	) ; done

ptest2: $(PROGS)
	set -e ; for i in $(PROGS); do ( \
	$(TEST_MPIRUN_2) ./$$i $(TEST_OUTDIR)/testfile.nc ; \
	) ; done

ptests: ptest2 ptest4
ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

