#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

LDADD  = ../common/libtestutils.la $(top_builddir)/src/libs/libpnetcdf.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

AM_CFLAGS =
if IS_BIGENDIAN
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   AM_CFLAGS += -DWORDS_BIGENDIAN
   # AM_FFLAGS  += $(FC_DEFINE)WORDS_BIGENDIAN
   # AM_FCFLAGS += $(FC_DEFINE)WORDS_BIGENDIAN
endif

check_PROGRAMS = bb_bsize \
                 bb_hints \
                 bb_many_reqs \
                 bb_nonblocking \
                 highdim \
				 varn

EXTRA_DIST = seq_runs.sh parallel_run.sh

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
			 $(TESTOUTDIR)/*.nc $(TESTOUTDIR)/*bb \
			 $(TESTOUTDIR)/*.meta $(TESTOUTDIR)/*data

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = export check_PROGRAMS="$(check_PROGRAMS)";
# AM_TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
# AM_TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";

if PNETCDF_DEBUG
   TESTS_ENVIRONMENT += export PNETCDF_DEBUG=1;
endif
if ENABLE_BURST_BUFFER
   TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER=1;
endif
if ENABLE_THREAD_SAFE
    TESTS_ENVIRONMENT += export ENABLE_THREAD_SAFE=1;
endif
if MIMIC_LUSTRE
   TESTS_ENVIRONMENT += export MIMIC_LUSTRE=1;
endif

TESTS = $(check_PROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/seq_runs.sh
SH_LOG_COMPILER =

# Some of these tests are designed to run on one processes,
# Run them on 4 processes to see if they can handle well

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 2 || exit 1

ptest6: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 6 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 6 || exit 1

ptests: ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

