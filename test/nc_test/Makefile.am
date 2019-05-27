#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .h .m4

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@ -lm

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

M4_HFILES = tests.h.m4

M4_SRCS   = test_get.m4 \
            test_put.m4 \
            test_iget.m4 \
            test_iput.m4 \
            test_read.m4 \
            test_write.m4

TESTPROGRAMS = t_nc \
               tst_misc \
               tst_norm \
               tst_small \
               tst_names \
               tst_atts3 \
               tst_atts

nc_test_SOURCES = nc_test.c error.c util.c
nodist_nc_test_SOURCES = tests.h $(M4_SRCS:.m4=.c)

BUILT_SOURCES = tests.h

$(M4_SRCS:.m4=.c): Makefile

M4FLAGS += -I${top_srcdir}/m4

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< > $@

tests.h: tests.h.m4 Makefile
	$(M4) $(M4FLAGS) $< > $@

check_PROGRAMS = $(TESTPROGRAMS) tst_nofill nc_test

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export PNETCDF_DEBUG="$(PNETCDF_DEBUG)";
TESTS_ENVIRONMENT += export TESTPROGRAMS="$(TESTPROGRAMS)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";
TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER="$(ENABLE_BURST_BUFFER)";
TESTS_ENVIRONMENT += export ENABLE_NETCDF4="$(ENABLE_NETCDF4)";

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = tests.h $(M4_SRCS:.m4=.c) \
             $(TESTOUTDIR)/test.nc \
             $(TESTOUTDIR)/scratch.nc \
             $(TESTOUTDIR)/tooth-fairy.nc \
             $(TESTOUTDIR)/tst_nofill.nc.fill \
             $(TESTOUTDIR)/tst_nofill.nc.nofill \
             $(TESTOUTDIR)/tst_atts3.nc.2 \
             $(TESTOUTDIR)/tst_atts3.bb.nc.2 \
             core core.* *.gcda *.gcno *.gcov gmon.out $(NC_FILES)

EXTRA_DIST = error.h $(M4_SRCS) $(M4_HFILES) README seq_runs.sh wrap_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# sequential runs only
ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

