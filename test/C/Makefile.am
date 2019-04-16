#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = pres_temp_4D_wr \
               pres_temp_4D_rd

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
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

TESTS = seq_runs.sh
TEST_EXTENSIONS = .sh
# LOG_COMPILER = $(srcdir)/wrap_runs.sh
# SH_LOG_COMPILER =

# Dependency: pres_temp_4D_rd reads the output file from pres_temp_4D_wr.
# To support parallel "make -jN check", we add dependencies below and automake
# states it currently works only for tests that end in one of the suffixes
# listed in TEST_EXTENSIONS. $(EXEEXT) are appended to work around.
pres_temp_4D_rd$(EXEEXT).log: pres_temp_4D_wr$(EXEEXT).log

EXTRA_DIST = seq_runs.sh parallel_run.sh

CLEANFILES = $(TESTOUTDIR)/pres_temp_4D.nc \
             $(TESTOUTDIR)/pres_temp_4D.nc4 \
             $(TESTOUTDIR)/pres_temp_4D.bb.nc \
             core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptests ptest4: $(TESTPROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

