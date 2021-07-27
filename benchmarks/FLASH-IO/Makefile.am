#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

ACLOCAL_AMFLAGS = -I m4

SUFFIXES = .o .F90

FC = $(MPIF90)

flash_benchmark_io_LDADD = -lpnetcdf

if BUILD_BENCHMARKS_IN_PNETCDF
   AM_FCFLAGS = $(FC_MODINC)${top_builddir}/src/binding/f90 $(FFREEFORMFLAG)
   AM_LDFLAGS = -L${top_builddir}/src/libs
   flash_benchmark_io_LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@
else
   AM_FCFLAGS = $(FC_MODINC)$(PNETCDF_INC)
   AM_LDFLAGS = $(PNETCDF_LIB)
endif

AM_FCFLAGS += -I$(srcdir)

AM_FCFLAGS += $(FC_DEFINE)N_DIM=3
AM_FCFLAGS += $(FC_DEFINE)MAXBLOCKS=100
AM_FCFLAGS += $(FC_DEFINE)IONMAX=13

check_PROGRAMS = flash_benchmark_io

noinst_PROGRAMS = $(check_PROGRAMS)

flash_benchmark_io_SOURCES = get_mfluid_property.F90 \
                             flash_release.F90 \
                             flash_benchmark_io.F90 \
                             checkpoint_ncmpi_parallel.F90 \
                             plotfile_ncmpi_parallel.F90

EXTRA_DIST = block_boundary_data.fh \
             common.fh \
             definitions.fh \
             numfluids.fh \
             physicaldata.fh \
             tree.fh \
             README \
             configure.ac \
             parallel_run.sh

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/flash_benchmark_io.ncmpi_chk_0000.nc \
             $(TESTOUTDIR)/flash_benchmark_io.ncmpi_plt_crn_0000.nc \
             $(TESTOUTDIR)/flash_benchmark_io.ncmpi_plt_cnt_0000.nc \
             $(TESTOUTDIR)/flash_benchmark_io.bb.ncmpi_chk_0000.nc \
             $(TESTOUTDIR)/flash_benchmark_io.bb.ncmpi_plt_crn_0000.nc \
             $(TESTOUTDIR)/flash_benchmark_io.bb.ncmpi_plt_cnt_0000.nc

TESTS =

TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export PNETCDF_DEBUG="$(PNETCDF_DEBUG)";
TESTS_ENVIRONMENT += export TESTPROGRAMS="$(TESTPROGRAMS)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";
TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER="$(ENABLE_BURST_BUFFER)";
if BUILD_BENCHMARKS_IN_PNETCDF
   TESTS_ENVIRONMENT += export BUILD_BENCHMARKS_IN_PNETCDF=1;
endif


ptest ptests ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

