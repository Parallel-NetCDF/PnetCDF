#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .f90

AM_DEFAULT_SOURCE_EXT = .f90

AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/binding/f90 \
             $(FC_MODINC)../common $(FFREEFORMFLAG)

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

TESTPROGRAMS = tst_f90 \
               f90tst_vars \
               tst_types2 \
               tst_f90_cdf5 \
               f90tst_vars2 \
               f90tst_vars3 \
               f90tst_vars4 \
               test_intent \
               test_attr_int64 \
               test_fill

PARALLEL_PROGS = f90tst_parallel \
                 f90tst_parallel2 \
                 f90tst_parallel3 \
                 f90tst_parallel4

check_PROGRAMS = $(TESTPROGRAMS) \
                 $(PARALLEL_PROGS) \
                 tst_io

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
TESTS_ENVIRONMENT += export PARALLEL_PROGS="$(PARALLEL_PROGS)";

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

EXTRA_DIST = seq_runs.sh wrap_runs.sh parallel_run.sh

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = *.$(FC_MODEXT) core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/tst_io1.nc $(TESTOUTDIR)/tst_io1.nc0 \
             $(TESTOUTDIR)/testfile.nc $(NC_FILES)

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptest4: $(PARALLEL_PROGS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2: $(PARALLEL_PROGS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 2 || exit 1

ptest8: $(PARALLEL_PROGS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 8 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 8 || exit 1

ptest10: $(PARALLEL_PROGS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 10 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 10 || exit 1

ptests: ptest2 ptest4 ptest8 ptest10
ptest6:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

