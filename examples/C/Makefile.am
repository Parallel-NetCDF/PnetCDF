#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c
AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS = -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la

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
                 req_all \
                 vard_mvars

# parallel runs only
# TESTS = $(check_PROGRAMS)

# Note: put_vara must be run immediately before get_vara

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(NC_FILES)

TESTMPIRUN3  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/3/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN8  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/8/g'`

ptest ptest4: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	for i in $(check_PROGRAMS); do { \
        if [ $$i = "get_vara" ] ; then \
           $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/put_vara.nc ; \
        else \
           $(TESTMPIRUN4) ./$$i -q $(TESTOUTDIR)/$$i.nc ; \
        fi ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done ; } ; done

ptest8: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	for i in $(check_PROGRAMS); do { \
        if [ $$i = "get_vara" ] ; then \
           $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/put_vara.nc ; \
        else \
           $(TESTMPIRUN8) ./$$i -q $(TESTOUTDIR)/$$i.nc ; \
        fi ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 8 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 8 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done ; } ; done

ptest3: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	for i in $(check_PROGRAMS) ; do { \
        if [ $$i = "get_vara" ] ; then \
           $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/put_vara.nc ; \
        else \
           $(TESTMPIRUN3) ./$$i -q $(TESTOUTDIR)/$$i.nc ; \
        fi ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:  C  parallel run on 3 processes --------------- $$i"; \
        else \
           echo "FAILED:  C  parallel run on 3 processes ------------- $$i"; \
	   exit 1; \
        fi ; } ; done ; } ; done

ptests: ptest3 ptest4 ptest8
ptest2 ptest6 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

