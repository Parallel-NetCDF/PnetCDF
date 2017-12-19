#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_FFLAGS    = -I$(top_builddir)/src/binding/f77
AM_FCFLAGS   = $(FC_MODINC)$(top_builddir)/src/binding/f90
AM_FCFLAGS  += $(FC_MODINC)$(srcdir)/../common
AM_FCFLAGS  += -I$(top_builddir)/src/binding/f77 
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la -lm

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS   = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = pnctest \
               csnap

if HAS_FORTRAN
TESTPROGRAMS += pnf_test pnctestf fixedform
pnf_test_SOURCES = pnf_test.f
pnctestf_SOURCES = pnctestf.f
fixedform_SOURCES = fixedform.f90
fixedform_FCFLAGS = $(FFIXEDFORMFLAG) $(FC_MODINC)$(top_builddir)/src/binding/f90 $(FC_MODINC)$(srcdir)/../common -I$(top_builddir)/src/binding/f77
if HAVE_F77_SUPPORT_FREEFORM
TESTPROGRAMS += freeform
freeform_SOURCES = freeform.f
freeform_FFLAGS = $(FFREEFORMFLAG) -I$(top_builddir)/src/binding/f77
endif
endif

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

# programs in this folder are just for testing compile error/warning messages,
# not for running
# TESTS = $(TESTPROGRAMS)

EXTRA_DIST = README

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

