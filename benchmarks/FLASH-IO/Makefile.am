#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

ACLOCAL_AMFLAGS = -I m4

SUFFIXES = .o .F90

FC          = $(MPIF90)

if BUILD_BENCHMARKS_IN_PNETCDF
PNETCDF_INC=${top_builddir}/src/binding/f90
PNETCDF_LIB="-L${top_builddir}/src/libs"
endif

AM_FCFLAGS  = $(FC_MODINC)$(PNETCDF_INC) -I$(srcdir) $(NAGf90FPPFLAGS)
AM_FCFLAGS += $(NAG_FCFLAGS)

AM_FCFLAGS += $(FC_DEFINE)N_DIM=3
AM_FCFLAGS += $(FC_DEFINE)MAXBLOCKS=100
AM_FCFLAGS += $(FC_DEFINE)IONMAX=13

if !BUILD_TESTSETS
noinst_PROGRAMS = flash_benchmark_io
else
RUN_QUIETLY = -q
endif
check_PROGRAMS = flash_benchmark_io

TESTS = $(check_PROGRAMS)

flash_benchmark_io_SOURCES = get_mfluid_property.F90 \
                             flash_release.F90 \
                             flash_benchmark_io.F90 \
                             checkpoint_ncmpi_parallel.F90 \
                             plotfile_ncmpi_parallel.F90

flash_benchmark_io_LDFLAGS = $(PNETCDF_LIB)
flash_benchmark_io_LDADD   = -lpnetcdf

EXTRA_DIST = block_boundary_data.fh \
             common.fh \
             definitions.fh \
             numfluids.fh \
             physicaldata.fh \
             tree.fh \
             README

CLEANFILES = *.nc core core.*

TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`

ptests ptest ptest4: $(check_PROGRAMS)
	$(TESTMPIRUN4) ./$(check_PROGRAMS) $(RUN_QUIETLY) $(TESTOUTDIR)/flash_io_test_ ; \
	if [ $$? = 0 ] ; then \
	    echo "PASS: F90 parallel run on 4 processes --------------- FLASH_IO"; \
	else \
	    echo "FAILED: F90 parallel run on 4 processes ------------- FLASH_IO"; \
	fi

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

