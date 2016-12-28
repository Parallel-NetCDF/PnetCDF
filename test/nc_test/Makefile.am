#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib -I$(srcdir)/../common
LDADD = $(top_builddir)/src/lib/libpnetcdf.la ../common/libtestutils.la -lm

M4_HFILES = tests.h.m4

M4_SRCS   = test_get.m4 \
            test_put.m4 \
            test_iget.m4 \
            test_iput.m4 \
            test_read.m4 \
            test_write.m4

TESTPROGRAMS = nc_test \
               t_nc \
               tst_misc \
               tst_norm \
               tst_small \
               tst_names \
               tst_atts3 \
               tst_atts \
               tst_nofill

nc_test_SOURCES = nc_test.c error.c util.c
nodist_nc_test_SOURCES = $(M4_SRCS:.m4=.c)

$(M4_SRCS:.m4=.c): Makefile

M4FLAGS += -I${top_srcdir}/m4

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

tests.h: tests.h.m4 Makefile
	$(M4) $(M4FLAGS) $< >$@

BUILT_SOURCES = tests.h $(M4_SRCS:.m4=.c)

check_PROGRAMS = $(TESTPROGRAMS)

AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

CLEANFILES = tests.h $(M4_SRCS:.m4=.c) \
             test.nc testfile.nc scratch.nc testfile.nc.2 \
             testfile.nc.nofill testfile.nc.fill core core.*

EXTRA_DIST = error.h $(M4_SRCS) $(M4_HFILES) README seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# sequential runs only
ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

