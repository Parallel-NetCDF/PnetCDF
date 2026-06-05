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
LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

check_PROGRAMS = collective_write \
                 nonblocking_write \
                 nonblocking_write_in_def \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 flexible_bottom \
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
                 transpose2D \
                 vard_int \
                 vard_bottom \
                 i_varn_int64 \
                 bput_varn_uint \
                 bput_varn_int64 \
                 fill_mode \
                 ghost_cell \
                 req_all \
                 vard_mvars \
                 time_var \
                 create_from_cdl

if INSTALL_EXAMPLES
   example_execbin_PROGRAMS = $(check_PROGRAMS)
   example_execbindir = $(exec_prefix)/pnetcdf_examples/C
   example_execbin_SCRIPTS = run_c_examples.sh
endif

if ENABLE_THREAD_SAFE
   check_PROGRAMS += pthread
endif

# parallel runs only
# TESTS = $(check_PROGRAMS)

# Note: put_vara must be run immediately before get_vara

TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export PNETCDF_DEBUG="$(PNETCDF_DEBUG)";
TESTS_ENVIRONMENT += export TESTPROGRAMS="$(TESTPROGRAMS)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";
TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER="$(ENABLE_BURST_BUFFER)";

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(NC_FILES) $(TESTOUTDIR)/pthread.nc.* $(TESTOUTDIR)/testfile.nc

EXTRA_DIST = run_c_examples.sh cdl_header.txt

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 4 "C " || exit 1

ptest8: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 8 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 8 "C " || exit 1

ptest3: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 3 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 3 "C " || exit 1

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

install-exec-hook:
	@if test "x$(example_execbindir)" != x ; then \
	cp -f $(srcdir)/cdl_header.txt $(DESTDIR)/$(example_execbindir) ; \
	$(SED_I) -e 's|check_PROGRAMS|$(check_PROGRAMS)|g ; s|TESTOUTDIR|$(TESTOUTDIR)|g ; s|TESTMPIRUN|$(TESTMPIRUN)|g ; s|SED_CMD|$(SED)|g ; s|ENABLE_BURST_BUFFER|$(ENABLE_BURST_BUFFER)|g ; s|ENABLE_NETCDF4|$(ENABLE_NETCDF4)|g' $(DESTDIR)/$(example_execbindir)/run_c_examples.sh ; \
	fi

