#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   # AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   # AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = header_consistency

check_PROGRAMS = $(TESTPROGRAMS)

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

if BUILD_DRIVER_BB
TESTBB = 1
TESTS_ENVIRONMENT += TESTBB=1 ; export TESTBB;
endif

# consistency check should be run on more than one process
# TESTS = seq_runs.sh

EXTRA_DIST = seq_runs.sh

CLEANFILES = $(TESTOUTDIR)/header_consistency.nc \
             core core.* *.gcda *.gcno *.gcov gmon.out

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`

ptest ptest4: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	if [ -n "$(TESTBB)" ]; then ( \
	export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable" ; \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	unset PNETCDF_HINTS ; \
	) ; fi ; \
	) ; done ; } ; done

ptest2: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	if [ -n "$(TESTBB)" ]; then ( \
	export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable" ; \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	unset PNETCDF_HINTS ; \
	) ; fi ; \
	) ; done ; } ; done

ptests: ptest2 ptest4
ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

