#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utils/ncmpidiff

LDADD  = ../common/libtestutils.la $(top_builddir)/src/utils/ncmpidiff/libncmpidiff_core.la
LDADD += $(top_builddir)/src/libs/libpnetcdf.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   # AM_FFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

check_PROGRAMS = test_inq_format \
                 cdf_type \
                 dim_cdf12 \
                 tst_open_cdf5 \
                 tst_corrupt

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = export check_PROGRAMS="$(check_PROGRAMS)";
# AM_TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
# AM_TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";
TESTS_ENVIRONMENT += export BAD_FILES="$(BAD_FILES)";

if PNETCDF_DEBUG
   TESTS_ENVIRONMENT += export PNETCDF_DEBUG=1;
endif
if ENABLE_BURST_BUFFER
   TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER=1;
endif
if ENABLE_THREAD_SAFE
    TESTS_ENVIRONMENT += export ENABLE_THREAD_SAFE=1;
endif
if MIMIC_LUSTRE
   TESTS_ENVIRONMENT += export MIMIC_LUSTRE=1;
endif

# During "make check", programs set in TESTS will be run by
# seq_runs.sh individually. Each produces its own log file.
TESTS = $(check_PROGRAMS) xfail_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/seq_runs.sh
SH_LOG_COMPILER =

BAD_FILES = bad_begin.nc5 \
            bad_xtype.nc1 bad_xtype.nc2 bad_xtype.nc5 \
            bad_ndims.nc1 bad_ndims.nc2 bad_ndims.nc5 \
            bad_dimid.nc1 bad_dimid.nc2 bad_dimid.nc5 \
            bad_nattrs.nc1 bad_nattrs.nc2 bad_nattrs.nc5

XFAIL_TESTS = xfail_runs.sh

CLEANFILES = $(TESTOUTDIR)/*.nc $(TESTOUTDIR)/*bb \
             core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = seq_runs.sh xfail_runs.sh parallel_run.sh \
             test_cdf.nc1 test_cdf.nc2 test_cdf.nc3 test_cdf.nc4 test_cdf.nc5 \
             $(BAD_FILES)

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

