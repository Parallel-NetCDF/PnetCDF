#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .F90 .c

AM_CPPFLAGS  = -I$(top_builddir)/src/include
AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_DEFAULT_SOURCE_EXT = .F90
AM_CFLAGS = -I$(top_srcdir)/src/binding/f77 -I$(top_builddir)/src/binding/f77
AM_FCFLAGS = -I$(top_builddir)/src/binding/f77 $(FC_MODINC)$(top_builddir)/src/binding/f90 $(FC_MODINC)../common
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la

if RELAX_COORD_BOUND
   AM_CPPFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
endif

M4SRCS      = test_get.m4 \
              test_put.m4  \
              test_iget.m4 \
              test_iput.m4

M4SRCS_F90  = $(M4SRCS:.m4=.F90)

F90_SRCS    = nf90_test.F90 \
              nf90_error.F90 \
              test_read.F90 \
              test_write.F90 \
              util.F90

TESTPROGRAMS = nf90_test

nf90_test_SOURCES = $(F90_SRCS) fortlib.c
nodist_nf90_test_SOURCES = $(M4SRCS_F90)

$(M4SRCS_F90): Makefile

M4FLAGS += -I${top_srcdir}/m4

.m4.F90:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

check_PROGRAMS = $(TESTPROGRAMS)

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

if ENABLE_NETCDF4
AM_CPPFLAGS += $(FC_DEFINE)ENABLE_NETCDF4
TESTNETCDF4 = 1
TESTS_ENVIRONMENT += TESTNETCDF4=1 ; export TESTNETCDF4;
endif

TESTS = seq_runs.sh
TEST_EXTENSIONS = .sh

if BUILD_DRIVER_BB
TESTBB = 1
TESTS_ENVIRONMENT += TESTBB=1 ; export TESTBB;
endif

CLEANFILES = $(M4SRCS_F90) \
             $(TESTOUTDIR)/scratch.nc \
             $(TESTOUTDIR)/test.nc \
             $(TESTOUTDIR)/tooth-fairy.nc \
             core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4SRCS) README seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# sequential runs only
ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

