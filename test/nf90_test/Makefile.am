#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .F90

AM_DEFAULT_SOURCE_EXT = .F90
AM_CFLAGS = -I$(top_srcdir)/src/libf -I$(top_builddir)/src/libf
AM_FCFLAGS = -I$(top_builddir)/src/libf $(FC_MODINC)$(top_builddir)/src/libf90 $(FC_MODINC)../common
LDADD = ${top_builddir}/src/lib/libpnetcdf.la ../common/libtestutils.la

if RELAX_COORD_BOUND
AM_FCFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
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

BUILT_SOURCES = $(M4SRCS_F90)

check_PROGRAMS = $(TESTPROGRAMS)

AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

CLEANFILES = $(M4SRCS_F90) scratch.nc test.nc core core.*

EXTRA_DIST = $(M4SRCS) tests.inc README seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# sequential runs only
ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

