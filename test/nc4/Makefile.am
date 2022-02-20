#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += @NETCDF4_INC@

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@ -lm

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS   = -DHAVE_DECL_MPI_OFFSET
endif

NCMPIGEN  = $(top_builddir)/src/utils/ncmpigen/ncmpigen
NCMPIDIFF = $(top_builddir)/src/utils/ncmpidiff/ncmpidiff

TESTPROGRAMS = simple_xy \
               put_get_all_kinds \
               interoperability_rd \
               interoperability_wr \
               noclobber \
               tst_zero_req \
               notsupport \
               tst_2_rec_dims \
               tst_get_put_size \
               pres_temp_4D \
               rd_compressed

if NETCDF_480_481
else
   TESTPROGRAMS += tst_rec_vars
endif

M4_SRCS  = put_get_all_kinds.m4 interoperability_rd.m4 interoperability_wr.m4

nodist_put_get_all_kinds_SOURCES = put_get_all_kinds.c
nodist_interoperability_rd_SOURCES = interoperability_rd.c
nodist_interoperability_wr_SOURCES = interoperability_wr.c

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

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

TESTS = $(TESTPROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = $(M4_SRCS:.m4=.c) \
             $(TESTOUTDIR)/put_get_all_kinds.nc.cdf4 \
             $(NC_FILES) testfile.nc

EXTRA_DIST = $(M4_SRCS) wrap_runs.sh parallel_run.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one process,
# Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# Run them on 2, 4, and 6 processes to see if they can handle well

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 2 || exit 1

ptest6: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 6 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 6 || exit 1

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

