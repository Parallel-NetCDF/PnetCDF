#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_FFLAGS    = -I$(top_builddir)/src/binding/f77
AM_FCFLAGS   = $(FC_MODINC)$(top_builddir)/src/binding/f90
AM_FCFLAGS  += $(FC_MODINC)$(srcdir)/../common
AM_FCFLAGS  += -I$(top_builddir)/src/binding/f77
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la -lm

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = test_bput \
               interleaved \
               i_varn_int64 \
               flexible_bput \
               wait_after_indep \
               req_all \
               i_varn_indef

M4_SRCS  = bput_varn.m4 \
           column_wise.m4

TESTPROGRAMS += $(M4_SRCS:.m4=)

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

M4FLAGS += -I${top_srcdir}/m4

if HAS_FORTRAN
   TESTPROGRAMS += mcoll_testf77 \
                   test_bputf77

   mcoll_testf77_SOURCES = mcoll_testf77.f
   test_bputf77_SOURCES = test_bputf77.f
if HAVE_MPI_MOD
   TESTPROGRAMS += mcoll_testf \
                   test_bputf
   mcoll_testf_SOURCES = mcoll_testf.f90
   test_bputf_SOURCES = test_bputf.f90
endif
endif

check_PROGRAMS = $(TESTPROGRAMS) mcoll_perf

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = seq_runs.sh

CLEANFILES = $(M4_SRCS:.m4=.c) *.nc core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_SRCS) seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one process,
# # Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# # Run them on 2, 4, and 6 processes to see if they can handle well

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN6  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/6/g'`

ptest ptest4: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/redef1.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptest2 : $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/redef1.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptest6 : $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/redef1.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN6) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

