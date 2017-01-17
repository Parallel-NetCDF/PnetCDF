#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS =
LDADD = $(top_builddir)/src/libs/libpnetcdf.la

if DECL_MPI_OFFSET
AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

check_PROGRAMS = aggregation \
                 write_block_read_column

# parallel runs only
# TESTS = $(check_PROGRAMS)

CLEANFILES = *.nc core core.*

ptest ptests ptest4: $(check_PROGRAMS)
	for i in $(check_PROGRAMS); do { \
	$(TESTMPIRUN4) ./$$i -q 10 $(TESTOUTDIR)/testfile.nc ; \
	if [ $$? = 0 ] ; then \
	    echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
	else \
	    echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
	    exit 1; \
	fi ; } ; done

ptest2 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

