#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

if IS_BIGENDIAN
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   AM_CFLAGS = -DWORDS_BIGENDIAN
   # AM_FFLAGS  += $(FC_DEFINE)WORDS_BIGENDIAN
   # AM_FCFLAGS += $(FC_DEFINE)WORDS_BIGENDIAN
endif

TESTPROGRAMS = large_files \
               large_var \
               large_dims_vars_attrs \
               high_dim_var \
               tst_cdf5_begin \
               large_coalesce

if HAS_FORTRAN
   TESTPROGRAMS += bigrecords
   bigrecords_SOURCES = bigrecords.f
   AM_FFLAGS = -I$(top_builddir)/src/binding/f77 $(FFIXEDFORMFLAG)
if HAVE_MPI_MOD
   TESTPROGRAMS += tst_flarge
   tst_flarge_SOURCES = tst_flarge.f90
   AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/binding/f90 \
                $(FC_MODINC)../common $(FFREEFORMFLAG)
endif
endif

EXTRA_DIST = wrap_runs.sh parallel_run.sh

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = $(NC_FILES) core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = export TESTPROGRAMS="$(TESTPROGRAMS)";
# AM_TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
# AM_TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
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

TESTS = $(check_PROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

# Some of these tests are designed to run on one processes,
# Run them on 4 processes to see if they can handle well

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptests: ptest4
ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(TESTPROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

