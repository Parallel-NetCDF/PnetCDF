#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .f90 .F90

AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS = -I$(top_builddir)/src/include

AM_FFLAGS = -I$(top_builddir)/src/binding/f77 $(FFIXEDFORMFLAG)
AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/binding/f90 $(FFREEFORMFLAG)
LDADD  = ${top_builddir}/src/libs/libpnetcdf.la
LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

check_PROGRAMS = pnetcdf-write-from-master \
                 pnetcdf-read-from-master \
                 pnetcdf-write-nfiles \
                 pnetcdf-read-nfiles \
                 pnetcdf-write-standard \
                 pnetcdf-read-standard \
                 pnetcdf-write-flexible \
                 pnetcdf-read-flexible \
                 pnetcdf-write-nb \
                 pnetcdf-read-nb \
                 pnetcdf-write-buffered \
                 pnetcdf-permute

if HAS_FORTRAN
check_PROGRAMS += pnetcdf-write-bufferedf77
pnetcdf_write_bufferedf77_SOURCES = pnetcdf-write-bufferedf77.f
pnetcdf_write_bufferedf77_LDADD = $(LDADD) ../F77/utils.o
if HAVE_MPI_MOD
check_PROGRAMS += pnetcdf-write-bufferedf
pnetcdf_write_bufferedf_SOURCES = pnetcdf-write-bufferedf.f90
pnetcdf_write_bufferedf_LDADD = $(LDADD) ../F77/utils.o
endif
endif

../F77/utils.o:
	set -e; cd ../F77 && $(MAKE) $(MFLAGS) utils.o

# parallel runs only
# TESTS = $(check_PROGRAMS)

TESTS_ENVIRONMENT  = export SED="$(SED)";
TESTS_ENVIRONMENT += export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";
TESTS_ENVIRONMENT += export TESTSEQRUN="$(TESTSEQRUN)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export PNETCDF_DEBUG="$(PNETCDF_DEBUG)";
TESTS_ENVIRONMENT += export TESTPROGRAMS="$(TESTPROGRAMS)";
TESTS_ENVIRONMENT += export check_PROGRAMS="$(check_PROGRAMS)";
TESTS_ENVIRONMENT += export ENABLE_BURST_BUFFER="$(ENABLE_BURST_BUFFER)";

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc) \
           $(check_PROGRAMS:%=$(TESTOUTDIR)/%.bb.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/pnetcdf-write-nfiles.*.nc \
             $(NC_FILES)

EXTRA_DIST = parallel_run.sh

ptest ptest4: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    $(subdir): Parallel testing on 4 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 4 || exit 1

ptests: ptest4
ptest2 ptest3 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

