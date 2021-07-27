#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .F

AM_DEFAULT_SOURCE_EXT = .F

AM_CPPFLAGS  = -I$(top_builddir)/src/include \
               -I$(top_srcdir)/src/include \
               -I$(top_builddir)/src/binding/f77 \
               -I$(top_srcdir)/src/binding/f77
AM_CFLAGS =
AM_FFLAGS = $(FFIXEDFORMFLAG)

LDADD  = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

M4FLAGS += -I${top_srcdir}/m4

if HAVE_F77_GNU_INT
   M4FFLAGS += -DHAVE_F77_GNU_INT
endif
if HAVE_F77_INT1
   M4FFLAGS += -DHAVE_F77_INT1
endif
if HAVE_F77_INT2
   M4FFLAGS += -DHAVE_F77_INT2
endif
if HAVE_F77_INT8
   M4FFLAGS += -DHAVE_F77_INT8
endif
if RELAX_COORD_BOUND
   AM_CFLAGS += -DRELAX_COORD_BOUND
   AM_FFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
endif
if ENABLE_NETCDF4
   AM_CFLAGS += -DENABLE_NETCDF4
   AM_FFLAGS += $(FC_DEFINE)ENABLE_NETCDF4
endif

M4SRCS = test_get.m4 \
         test_put.m4  \
         test_iget.m4 \
         test_iput.m4

.m4.F:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

F_SRCS = nf_test.F \
         nf_error.F \
         test_read.F \
         test_write.F \
         util.F

HFILES = tests.inc.in

TESTPROGRAMS = nf_test
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
TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER="$(ENABLE_BURST_BUFFER)";
TESTS_ENVIRONMENT += export ENABLE_NETCDF4="$(ENABLE_NETCDF4)";

TESTS = seq_runs.sh
TEST_EXTENSIONS = .sh

nf_test_SOURCES = $(F_SRCS) fortlib.c
nodist_nf_test_SOURCES = $(M4SRCS:.m4=.F)

$(M4SRCS:.m4=.F): Makefile

CLEANFILES = $(M4SRCS:.m4=.F) \
             $(TESTOUTDIR)/scratch.nc \
             $(TESTOUTDIR)/test.nc \
             $(TESTOUTDIR)/tooth-fairy.nc \
             core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4SRCS) $(HFILES) README seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# sequential runs only
ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

