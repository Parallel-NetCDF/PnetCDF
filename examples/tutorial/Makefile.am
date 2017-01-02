#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib
AM_FFLAGS = -I$(top_builddir)/src/libf
AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/libf90
LDADD = ${top_builddir}/src/lib/libpnetcdf.la

if DECL_MPI_OFFSET
AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
endif

check_PROGRAMS = pnetcdf-write-from-master \
                 pnetcdf-read-from-master \
                 pnetcdf-write-nfiles \
                 pnetcdf-read-nfiles \
                 pnetcdf-write-standard \
                 pnetcdf-read-standard \
                 pnetcdf-write-flexible \
                 pnetcdf-read-flexible \
                 pnetcdf-write-nb \
                 pnetcdf-read-nb \
                 pnetcdf-write-buffered \
                 pnetcdf-permute

if HAS_FORTRAN
check_PROGRAMS += pnetcdf-write-bufferedf77
pnetcdf_write_bufferedf77_SOURCES = pnetcdf-write-bufferedf77.f
if HAVE_MPI_MOD
check_PROGRAMS += pnetcdf-write-bufferedf
pnetcdf_write_bufferedf_SOURCES = pnetcdf-write-bufferedf.f90
endif
endif

# parallel runs only
# TESTS = $(check_PROGRAMS)

CLEANFILES = *.nc core core.*

TESTMPIRUN4  = `echo $(TESTMPIRUN) | sed -e 's/NP/4/g'`

ptest ptest4: $(check_PROGRAMS)
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:     parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED:     parallel run on 4 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done

ptests: ptest4
ptest2 ptest3 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests: $(check_PROGRAMS)
.PHONY: tests ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

