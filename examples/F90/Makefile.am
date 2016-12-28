#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@


SUFFIXES = .o .f .F90 .f90

AM_DEFAULT_SOURCE_EXT = .f90
AM_FCFLAGS = ${FC_MODINC}$(top_builddir)/src/libf90
LDADD = ${top_builddir}/src/lib/libpnetcdf.la

TESTPROGRAMS = nonblocking_write \
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

nonblocking_write_SOURCES = nonblocking_write.f90 utils.F90
column_wise_SOURCES = column_wise.f90 utils.F90
block_cyclic_SOURCES = block_cyclic.f90 utils.F90
flexible_api_SOURCES = flexible_api.f90 utils.F90
get_info_SOURCES = get_info.f90 utils.F90
hints_SOURCES = hints.f90 utils.F90
put_var_SOURCES =  put_var.f90 utils.F90
put_varn_real_SOURCES = put_varn_real.f90 utils.F90
put_varn_int_SOURCES = put_varn_int.f90 utils.F90
transpose_SOURCES = transpose.f90 utils.F90
vard_int_SOURCES = vard_int.f90 utils.F90
fill_mode_SOURCES = fill_mode.f90 utils.F90

check_PROGRAMS = $(TESTPROGRAMS)

# parallel runs only
# TESTS = $(TESTPROGRAMS)

CLEANFILES = *.nc core core.*

TESTMPIRUN4  = `echo $(TESTMPIRUN) | sed -e 's/NP/4/g'`
TESTMPIRUN8  = `echo $(TESTMPIRUN) | sed -e 's/NP/8/g'`
TESTMPIRUN3  = `echo $(TESTMPIRUN) | sed -e 's/NP/3/g'`

ptest ptest4: $(TESTPROGRAMS)
	@for i in $(TESTPROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 4 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptest8: $(TESTPROGRAMS)
	@for i in $(TESTPROGRAMS); do { \
        $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 8 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 8 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptest3: $(TESTPROGRAMS)
	@for i in $(TESTPROGRAMS); do { \
        $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F90 parallel run on 3 processes --------------- $$i"; \
        else \
           echo "FAILED: F90 parallel run on 3 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

