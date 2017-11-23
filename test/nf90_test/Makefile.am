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
   AM_FCFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
endif
if DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
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
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

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

