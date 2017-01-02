#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib -I$(srcdir)/../common
LDADD = $(top_builddir)/src/lib/libpnetcdf.la ../common/libtestutils.la

if IS_BIGENDIAN
   AM_CPPFLAGS += -DWORDS_BIGENDIAN
endif

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = large_files large_var

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

EXTRA_DIST = seq_runs.sh

CLEANFILES = *.nc core core.*

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one processes,
# Run them on 4 processes to see if they can handle well

ptest ptest4: $(TESTPROGRAMS)
	$(RM) -f $(TEST_OUTDIR)/testfile.nc $(TEST_OUTDIR)/testfile.nc
	set -e ; for i in $(TESTPROGRAMS) ; do ( \
	$(TEST_MPIRUN_4) ./$$i $(TEST_OUTDIR)/testfile.nc ; \
	) ; done

ptest2: $(TESTPROGRAMS)
	$(RM) -f $(TEST_OUTDIR)/testfile.nc $(TEST_OUTDIR)/testfile.nc
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TEST_MPIRUN_2) ./$$i $(TEST_OUTDIR)/testfile.nc ; \
	) ; done

ptest6: $(TESTPROGRAMS)
	$(RM) -f $(TEST_OUTDIR)/testfile.nc $(TEST_OUTDIR)/testfile.nc
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TEST_MPIRUN_6) ./$$i $(TEST_OUTDIR)/testfile.nc ; \
	) ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

