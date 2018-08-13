#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   # AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   # AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = header_consistency

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = export TESTPROGRAMS="$(TESTPROGRAMS)";
# AM_TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
# AM_TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT  = export TESTPROGRAMS="$(TESTPROGRAMS)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";

if BUILD_DRIVER_BB
   TESTS_ENVIRONMENT += export TESTBB=1;
endif

# consistency check should be run on more than one process
# TESTS = seq_runs.sh

EXTRA_DIST = seq_runs.sh parallel_run.sh

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = $(NC_FILES) core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export SED="$(SED)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	./parallel_run.sh 4

ptest2: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	./parallel_run.sh 2

ptests: ptest2 ptest4
ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

