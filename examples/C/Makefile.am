#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib
LDADD = $(top_builddir)/src/lib/libpnetcdf.la

check_PROGRAMS = collective_write \
                 nonblocking_write \
                 nonblocking_write_in_def \
                 column_wise \
                 block_cyclic \
                 flexible_api \
                 get_info \
                 hints \
                 mput \
                 put_varn_float \
                 put_varn_int \
                 create_open \
                 global_attributes \
                 put_vara \
                 get_vara \
                 transpose \
                 vard_int \
                 i_varn_int64 \
                 bput_varn_uint \
                 bput_varn_int64 \
                 fill_mode \
                 ghost_cell \
                 req_all

# parallel runs only
# TESTS = $(check_PROGRAMS)

# Note: put_vara must be run immediately before get_vara

CLEANFILES = *.nc core core.*

ptest ptest4: $(check_PROGRAMS)
	for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done

ptest8: $(check_PROGRAMS)
	for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 8 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 8 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done

ptest3: $(check_PROGRAMS)
	for i in $(check_PROGRAMS) ; do { \
        $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 3 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 3 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

