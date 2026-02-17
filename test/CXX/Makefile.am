#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .cpp

AM_DEFAULT_SOURCE_EXT = .cpp

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/binding/cxx
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include

LDADD  = ../common/libtestutils.la $(top_builddir)/src/libs/libpnetcdf.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

check_PROGRAMS = nctst \
                 test_classic

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";

if PNETCDF_DEBUG
   TESTS_ENVIRONMENT += export PNETCDF_DEBUG=1;
endif
if ENABLE_BURST_BUFFER
   TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER=1;
endif
if ENABLE_NETCDF4
   TESTS_ENVIRONMENT += export ENABLE_NETCDF4=1;
endif
if ENABLE_THREAD_SAFE
    TESTS_ENVIRONMENT += export ENABLE_THREAD_SAFE=1;
endif
if MIMIC_LUSTRE
   TESTS_ENVIRONMENT += export MIMIC_LUSTRE=1;
endif

TESTS = $(check_PROGRAMS)
LOG_COMPILER = $(srcdir)/seq_runs.sh

EXTRA_DIST = seq_runs.sh parallel_run.sh

CLEANFILES = $(TESTOUTDIR)/*.nc $(TESTOUTDIR)/*bb \
             core core.* *.gcda *.gcno *.gcov gmon.out a.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests


ptest ptests ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

