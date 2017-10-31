#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
LDADD += ../common/libtestutils.la

if ALL_STATIC
   AM_LDFLAGS = -all-static -static
endif

AM_CFLAGS =
if IS_BIGENDIAN
   AM_CFLAGS += -DWORDS_BIGENDIAN
endif

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS += -DHAVE_DECL_MPI_OFFSET
   # AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   # AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = large_files \
               large_var \
               large_dims_vars_attrs \
               high_dim_var \
               tst_cdf5_begin

check_PROGRAMS = $(TESTPROGRAMS)

EXTRA_DIST = seq_runs.sh

CLEANFILES = *.nc core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

TESTS =

if RUN_LARGE_FILE_TEST
TESTS += seq_runs.sh

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

# Some of these tests are designed to run on one processes,
# Run them on 4 processes to see if they can handle well

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN6  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/6/g'`

ptest ptest4: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/testfile.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS) ; do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptest2: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/testfile.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptest6: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	$(RM) -f $(TESTOUTDIR)/testfile.nc $(TESTOUTDIR)/testfile.nc ; \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN6) ./$$i $(TESTOUTDIR)/testfile.nc ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest6
else
ptest ptest2 ptest4 ptest6 ptests:
endif
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

