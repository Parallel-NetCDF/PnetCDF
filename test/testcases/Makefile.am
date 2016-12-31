#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib -I$(srcdir)/../common
AM_FFLAGS = -I$(top_builddir)/src/libf
AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/libf90 $(FC_MODINC)$(srcdir)/../common -I$(top_builddir)/src/libf 
LDADD = ${top_builddir}/src/lib/libpnetcdf.la ../common/libtestutils.la -lm

if SIZEOF_MPI_AINT_IS_4
AM_FFLAGS += $(FC_DEFINE)SIZEOF_MPI_AINT_IS_4
endif

NCMPIGEN  = $(top_builddir)/src/utils/ncmpigen/ncmpigen
NCMPIDIFF = $(top_builddir)/src/utils/ncmpidiff/ncmpidiff

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
               varn_contig \
               ivarn \
               check_striping \
               add_var \
               buftype_free \
               last_large_var \
               check_type \
               test_erange \
               scalar

M4_SRCS  = erange_fill.m4

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

if ENABLE_ERANGE_FILL
   M4FLAGS += -DERANGE_FILL
   TESTPROGRAMS += erange_fill
   BUILT_SOURCES = erange_fill.c
endif

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
if BUILD_LARGE_FILE_TEST
   TESTPROGRAMS += bigrecords
   bigrecords_SOURCES = bigrecords.f
endif
endif

check_PROGRAMS = $(TESTPROGRAMS) redef1

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

CLEANFILES = $(M4_SRCS:.m4=.c) *.nc core core.*

EXTRA_DIST = $(M4_SRCS) seq_runs.sh redef-good.ncdump

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one process,
# Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# Run them on 2, 4, and 6 processes to see if they can handle well

ptest ptest4: $(TESTPROGRAMS)
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done

ptest2: $(TESTPROGRAMS)
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done

ptest6: $(TESTPROGRAMS)
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN_6) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

