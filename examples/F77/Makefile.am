#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .f .F .F90

AM_DEFAULT_SOURCE_EXT = .f
AM_FFLAGS = -I$(top_builddir)/src/libf
LDADD = ${top_builddir}/src/lib/libpnetcdf.la

check_PROGRAMS = nonblocking_write \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 get_info \
                 hints \
                 put_vara \
                 put_varn_real \
                 put_varn_int \
                 transpose \
                 i_varn_real \
                 bput_varn_int8 \
                 fill_mode \
                 vard_int

nonblocking_write_SOURCES = nonblocking_write.f utils.F90
column_wise_SOURCES = column_wise.f utils.F90
block_cyclic_SOURCES = block_cyclic.f utils.F90
flexible_api_SOURCES = flexible_api.f utils.F90
get_info_SOURCES = get_info.f utils.F90
hints_SOURCES = hints.f utils.F90
put_vara_SOURCES =  put_vara.f utils.F90
put_varn_real_SOURCES = put_varn_real.f utils.F90
put_varn_int_SOURCES = put_varn_int.f utils.F90
transpose_SOURCES = transpose.f utils.F90
i_varn_real_SOURCES = i_varn_real.f utils.F90
bput_varn_int8_SOURCES = bput_varn_int8.f utils.F90
fill_mode_SOURCES = fill_mode.f utils.F90
vard_int_SOURCES = vard_int.F utils.F90

# parallel runs only
# TESTS = $(check_PROGRAMS)

CLEANFILES = *.nc core core.*

ptest ptest4: $(check_PROGRAMS)
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F77 parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED: F77 parallel run on 4 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptest8: $(check_PROGRAMS)
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F77 parallel run on 8 processes --------------- $$i"; \
        else \
           echo "FAILED: F77 parallel run on 8 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptest3: $(check_PROGRAMS)
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS: F77 parallel run on 3 processes --------------- $$i"; \
        else \
           echo "FAILED: F77 parallel run on 3 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

