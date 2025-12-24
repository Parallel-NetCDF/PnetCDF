#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .f .F90

AM_DEFAULT_SOURCE_EXT = .f
AM_FFLAGS = -I$(top_builddir)/src/binding/f77 $(FFIXEDFORMFLAG)
LDADD  = ${top_builddir}/src/libs/libpnetcdf.la utils.o
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

# suppress warning messages when using NAG Fortran compiler
if NAGFORT
  # suppresses warning messages about unused external procedures
  AM_FFLAGS += -w=uep

  # suppresses extension warnings for obsolete but common extensions to Fortran
  # 77 - these are TAB format, byte-length specifiers, Hollerith constants
  # and D lines;
  AM_FFLAGS += -w=x77

  # suppresses warning messages about variables set but never referenced;
  AM_FFLAGS += -w=unreffed

  # suppresses warning messages about unused PARAMETERs;
  AM_FFLAGS += -w=uparam
endif

if SIZEOF_MPI_AINT_IS_4
   AM_FFLAGS += $(FC_DEFINE)SIZEOF_MPI_AINT_IS_4
endif

check_PROGRAMS = nonblocking_write \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 get_info \
                 hints \
                 put_vara \
                 get_vara \
                 put_varn_real \
                 put_varn_int \
                 transpose \
                 i_varn_real \
                 bput_varn_int8 \
                 fill_mode \
                 vard_int \
                 time_var

if INSTALL_EXAMPLES
   example_execbin_PROGRAMS = $(check_PROGRAMS)
   example_execbindir = $(exec_prefix)/pnetcdf_examples/F77
   example_execbin_SCRIPTS = run_f77_examples.sh
endif

EXTRA_DIST = utils.F90 run_f77_examples.sh

utils.o: $(srcdir)/utils.F90
	$(FC) $(FC_DEFS) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ $<

# parallel runs only
# TESTS = $(check_PROGRAMS)

TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(FSTYPE_PREFIX)$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export TESTPROGRAMS="$(TESTPROGRAMS)";
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

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(NC_FILES)

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 4 "F77" || exit 1

ptest8: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 8 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 8 "F77" || exit 1

ptest3: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 3 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/../parallel_run.sh 3 "F77" || exit 1

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

install-exec-hook:
	@if test "x$(example_execbindir)" != x ; then \
	$(SED_I) -e 's|check_PROGRAMS|$(check_PROGRAMS)|g ; s|TESTOUTDIR|$(TESTOUTDIR)|g ; s|TESTMPIRUN|$(TESTMPIRUN)|g ; s|SED_CMD|$(SED)|g ; s|ENABLE_BURST_BUFFER|$(ENABLE_BURST_BUFFER)|g ; s|ENABLE_NETCDF4|$(ENABLE_NETCDF4)|g' $(DESTDIR)/$(example_execbindir)/run_f77_examples.sh ; \
	fi

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

