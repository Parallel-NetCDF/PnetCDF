#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la

if DECL_MPI_OFFSET
   AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = pres_temp_4D_wr \
               pres_temp_4D_rd

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
TESTS_ENVIRONMENT  = TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = $(check_PROGRAMS)
LOG_COMPILER = seq_runs.sh

# Dependency: pres_temp_4D_rd reads the output file from pres_temp_4D_wr.
# To support parallel "make -jN check", we add dependencies below and automake
# states it currently works only for tests that end in one of the suffixes
# listed in TEST_EXTENSIONS. $(EXEEXT) are appended to work around.
pres_temp_4D_rd$(EXEEXT).log: pres_temp_4D_wr$(EXEEXT).log

EXTRA_DIST = seq_runs.sh

CLEANFILES = $(TESTOUTDIR)/pres_temp_4D.nc \
             core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`

ptest ptests ptest4: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/pres_temp_4D.nc ; \
	) ; done ; } ; done

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

