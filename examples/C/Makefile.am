#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c
AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS = -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la

check_PROGRAMS = collective_write \
                 nonblocking_write \
                 nonblocking_write_in_def \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 get_info \
                 hints \
                 mput \
                 put_varn_float \
                 put_varn_int \
                 create_open \
                 global_attributes \
                 put_vara \
                 get_vara \
                 transpose \
                 vard_int \
                 i_varn_int64 \
                 bput_varn_uint \
                 bput_varn_int64 \
                 fill_mode \
                 ghost_cell \
                 req_all \
                 vard_mvars

if ENABLE_THREAD_SAFE
   check_PROGRAMS += pthread
endif

# parallel runs only
# TESTS = $(check_PROGRAMS)

# Note: put_vara must be run immediately before get_vara

TESTS_ENVIRONMENT =
if BUILD_DRIVER_BB
   TESTS_ENVIRONMENT += export TESTBB=1;
endif

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(NC_FILES) $(TESTOUTDIR)/pthread.nc.*

EXTRA_DIST = parallel_run.sh

TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export SED="$(SED)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	./parallel_run.sh 4

ptest8: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 8 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	./parallel_run.sh 8

ptest3: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 3 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	./parallel_run.sh 3

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

