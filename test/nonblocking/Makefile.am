#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_srcdir)/src/include \
              -I$(srcdir)/../common \
              -I$(top_builddir)/src/include \
              -I$(top_builddir)/src/binding/f77

AM_FCFLAGS  = $(FC_MODINC)$(top_builddir)/src/binding/f90 \
              $(FC_MODINC)$(srcdir)/../common $(FFREEFORMFLAG)
AM_FFLAGS   = -I$(top_builddir)/src/binding/f77 $(FFIXEDFORMFLAG)

LDADD  = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@ -lm

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = test_bput \
               interleaved \
               i_varn_int64 \
               flexible_bput \
               wait_after_indep \
               req_all \
               i_varn_indef \
               large_num_reqs

M4_SRCS      = bput_varn.m4 \
               column_wise.m4

TESTPROGRAMS += $(M4_SRCS:.m4=)

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

M4FLAGS += -I${top_srcdir}/m4

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# # and 'make install'. This means you cannot build a specific target (e.g.,
# # 'make foo') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

nodist_bput_varn_SOURCES = bput_varn.c
nodist_column_wise_SOURCES = column_wise.c

if HAS_FORTRAN
   TESTPROGRAMS += mcoll_testf77 \
                   test_bputf77

   mcoll_testf77_SOURCES = mcoll_testf77.f
   test_bputf77_SOURCES = test_bputf77.f
if HAVE_MPI_MOD
   TESTPROGRAMS += mcoll_testf \
                   test_bputf
   mcoll_testf_SOURCES = mcoll_testf.f90
   test_bputf_SOURCES = test_bputf.f90
endif
endif

check_PROGRAMS = $(TESTPROGRAMS) mcoll_perf

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

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/testfile*.nc $(NC_FILES) \
             $(TESTOUTDIR)/mcoll_perf.nc.* $(TESTOUTDIR)/mcoll_perf.bb.nc.*

EXTRA_DIST = $(M4_SRCS) seq_runs.sh wrap_runs.sh parallel_run.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one process,
# # Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# # Run them on 2, 4, and 6 processes to see if they can handle well

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2 : $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 2 || exit 1

ptest6 : $(check_PROGRAMS)
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

