#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS = -I$(top_builddir)/src/include
LDADD = $(top_builddir)/src/libs/libpnetcdf.la

check_PROGRAMS = aggregation \
                 write_block_read_column

# parallel runs only
# TESTS = $(check_PROGRAMS)

NC_FILES = $(check_PROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(NC_FILES)

TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`

if BUILD_DRIVER_BB
TESTBB = 1
TESTS_ENVIRONMENT = TESTBB=1 ; export TESTBB;
endif

ptest ptests ptest4: $(check_PROGRAMS)
	for i in $(check_PROGRAMS); do { \
	$(TESTMPIRUN4) ./$$i -q 10 $(TESTOUTDIR)/$$i.nc ; \
	if [ $$? = 0 ] ; then \
	    echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
	else \
	    echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
	    exit 1; \
	fi ; \ 
	if [ -n "$(TESTBB)" ]; then ( \
        export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable" ; \
        $(TESTMPIRUN4) ./$$i -q 10 $(TESTOUTDIR)/$$i.nc ; \
        unset PNETCDF_HINTS ; \
		if [ $$? = 0 ] ; then \
			echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
		else \
			echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
			exit 1; \
		fi ; \ 
	) ; fi ; \
	} ; done

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

