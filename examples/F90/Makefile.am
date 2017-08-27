#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@


SUFFIXES = .o .f .F90 .f90

AM_DEFAULT_SOURCE_EXT = .f90
AM_FCFLAGS = ${FC_MODINC}$(top_builddir)/src/binding/f90
LDADD = ${top_builddir}/src/libs/libpnetcdf.la utils.o

check_PROGRAMS = nonblocking_write \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 get_info \
                 hints \
                 put_var \
                 put_varn_real \
                 put_varn_int \
                 transpose \
                 fill_mode \
                 vard_int

EXTRA_DIST = utils.F90

utils.o: utils.F90
	$(FC) $(FC_DEFS) -c -o $@ $(FCFLAGS) $<

# parallel runs only
# TESTS = $(check_PROGRAMS)

CLEANFILES = *.nc core core.* *.gcda *.gcno *.gcov gmon.out

TESTMPIRUN3  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/3/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN8  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/8/g'`

ptest ptest4: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 4 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done ; } ; done

ptest8: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 8 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 8 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done ; } ; done

ptest3: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 3 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 3 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done ; } ; done

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

