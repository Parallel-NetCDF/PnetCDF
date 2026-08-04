#
# Copyright (C) 2026, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS =

LDADD =

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

check_PROGRAMS =

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export top_builddir="$(top_builddir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";

TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER=@ENABLE_BURST_BUFFER@;
TESTS_ENVIRONMENT += export ENABLE_NETCDF4=@ENABLE_NETCDF4@;
TESTS_ENVIRONMENT += export ENABLE_THREAD_SAFE=@ENABLE_THREAD_SAFE@;
TESTS_ENVIRONMENT += export MIMIC_LUSTRE=@MIMIC_LUSTRE@;

# build check targets but not invoke by not setting TESTS
# TESTS = $(check_PROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/seq_runs.sh
SH_LOG_COMPILER =

check_SCRIPTS = seq_runs.sh parallel_run.sh

EXTRA_DIST = seq_runs.sh parallel_run.sh

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

