#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .f90

AM_DEFAULT_SOURCE_EXT = .f90
AM_FCFLAGS = -I$(top_builddir)/src/binding/f77 $(FC_MODINC)$(top_builddir)/src/binding/f90 $(FC_MODINC)../common
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la

TESTPROGRAMS = tst_f90 \
               f90tst_vars \
               tst_types2 \
               tst_f90_cdf5 \
               f90tst_vars2 \
               f90tst_vars3 \
               f90tst_vars4 \
               test_intent \
               test_attr_int64

PARALLEL_PROGS = f90tst_parallel \
                 f90tst_parallel2 \
                 f90tst_parallel3 \
                 f90tst_parallel4

check_PROGRAMS = $(TESTPROGRAMS) \
                 $(PARALLEL_PROGS) \
                 tst_io \
                 tst_flarge

TESTS_ENVIRONMENT =

if RUN_LARGE_FILE_TEST
   TESTPROGRAMS += tst_flarge
   TESTS_ENVIRONMENT += RUN_LARGE_FILE_TEST="yes" ; export RUN_LARGE_FILE_TEST;
endif

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT += TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

EXTRA_DIST = seq_runs.sh wrap_runs.sh

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = *.$(FC_MODEXT) core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/tst_io1.nc $(NC_FILES)

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN8  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/8/g'`
TESTMPIRUN10 = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/10/g'`


ptest ptest4: $(PARALLEL_PROGS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(PARALLEL_PROGS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest2: $(PARALLEL_PROGS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(PARALLEL_PROGS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest8: $(PARALLEL_PROGS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(PARALLEL_PROGS); do ( \
	$(TESTMPIRUN8) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest10: $(PARALLEL_PROGS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(PARALLEL_PROGS); do ( \
	$(TESTMPIRUN10) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest8 ptest10
ptest6:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

