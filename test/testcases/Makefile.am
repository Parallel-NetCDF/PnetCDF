#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_srcdir)/src/include \
              -I$(srcdir)/../common \
              -I$(top_builddir)/src/include \
              -I$(top_builddir)/src/binding/f77

AM_FCFLAGS  = $(FC_MODINC)$(top_builddir)/src/binding/f90 \
              $(FC_MODINC)$(srcdir)/../common $(FFREEFORMFLAG)
AM_FFLAGS   = -I$(top_builddir)/src/binding/f77 $(FFIXEDFORMFLAG)

LDADD       = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la
LDADD      += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@ -lm

if SIZEOF_MPI_AINT_IS_4
   AM_FFLAGS += $(FC_DEFINE)SIZEOF_MPI_AINT_IS_4
endif

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS   = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = ncmpi_vars_null_stride \
               vectors \
               collective_error \
               test_varm \
               alignment_test \
               flexible \
               flexible2 \
               flexible_varm \
               nonblocking \
               noclobber \
               record \
               inq_num_vars \
               varn_int \
               modes \
               one_record \
               inq_recsize \
               test_vard \
               test_vard_rec \
               test_vard_multiple \
               varn_contig \
               ivarn \
               check_striping \
               add_var \
               buftype_free \
               last_large_var \
               check_type \
               test_erange \
               scalar \
               null_args \
               tst_dimsizes \
               mix_collectives \
               large_var_cdf5 \
               tst_max_var_dims \
               tst_info \
               tst_vars_fill \
               tst_def_var_fill \
               test_fillvalue \
               error_precedence \
               tst_free_comm

M4_SRCS      = put_all_kinds.m4 \
               erange_fill.m4 \
               tst_vars_fill.m4 \
               iput_all_kinds.m4 \
               error_precedence.m4 \
               null_args.m4

nodist_put_all_kinds_SOURCES = put_all_kinds.c
nodist_erange_fill_SOURCES = erange_fill.c
nodist_tst_vars_fill_SOURCES = tst_vars_fill.c
nodist_iput_all_kinds_SOURCES = iput_all_kinds.c
nodist_error_precedence_SOURCES = error_precedence.c
nodist_null_args_SOURCES = null_args.c

if ENABLE_ERANGE_FILL
   M4FLAGS += -DERANGE_FILL
   TESTPROGRAMS += erange_fill
endif

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

if HAS_FORTRAN
   TESTPROGRAMS += varn_intf \
                   attrf \
                   buftype_freef \
                   put_parameter \
                   test_vardf

       varn_intf_SOURCES = varn_intf.f
           attrf_SOURCES = attrf.f
   buftype_freef_SOURCES = buftype_freef.f
   put_parameter_SOURCES = put_parameter.f
      test_vardf_SOURCES = test_vardf.F
if HAVE_MPI_MOD
   TESTPROGRAMS += inq_num_varsf \
                   inq_recsizef \
                   test_vardf90 \
                   varn_real

   inq_num_varsf_SOURCES = inq_num_varsf.f90
    inq_recsizef_SOURCES = inq_recsizef.f90
    test_vardf90_SOURCES = test_vardf90.f90
       varn_real_SOURCES = varn_real.f90
endif
endif

# all programs in TESTPROGRAMS will be run by wrap_runs.sh
# others in check_PROGRAMS but not in TESTPROGRAMS will be run by seq_runs.sh
# Those are the ones need special treatment
check_PROGRAMS = $(TESTPROGRAMS) \
                 put_all_kinds \
                 redef1 \
                 iput_all_kinds \
                 tst_version

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
TESTS_ENVIRONMENT += export ENABLE_THREAD_SAFE="$(ENABLE_THREAD_SAFE)";

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc1) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc2) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc3) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc4) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc5) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc1) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc2) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc3) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc4) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc5)

if ENABLE_THREAD_SAFE
   check_PROGRAMS += tst_pthread
endif

BURST_BUFFER_FILES = $(NC_FILES:=_*.meta) $(NC_FILES:=_*.data)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/redef1.nc $(TESTOUTDIR)/redef1.bb.nc \
             $(TESTOUTDIR)/redef2.nc \
             $(TESTOUTDIR)/tst_pthread.nc.* $(TESTOUTDIR)/testfile.nc* \
             $(NC_FILES)

EXTRA_DIST = $(M4_SRCS) seq_runs.sh redef-good.ncdump \
             wrap_runs.sh parallel_run.sh

# Some of these tests are designed to run on one process,
# Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# Run them on 2, 4, and 6 processes to see if they can handle well

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 2 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 2 || exit 1

ptest6: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 6 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 6 || exit 1

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

