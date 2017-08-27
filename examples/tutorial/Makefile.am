#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .f90 .F90

AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS = -I$(top_builddir)/src/include

AM_FFLAGS = -I$(top_builddir)/src/binding/f77
AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/src/binding/f90
LDADD = ${top_builddir}/src/libs/libpnetcdf.la

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
pnetcdf_write_bufferedf77_LDADD = $(LDADD) ../F77/utils.o
if HAVE_MPI_MOD
check_PROGRAMS += pnetcdf-write-bufferedf
pnetcdf_write_bufferedf_SOURCES = pnetcdf-write-bufferedf.f90
pnetcdf_write_bufferedf_LDADD = $(LDADD) ../F77/utils.o
endif
endif

../F77/utils.o:
	set -e; cd ../F77 && $(MAKE) $(MFLAGS) utils.o

# parallel runs only
# TESTS = $(check_PROGRAMS)

CLEANFILES = *.nc core core.* *.gcda *.gcno *.gcov gmon.out

TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`

ptest ptest4: $(check_PROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	@for i in $(check_PROGRAMS); do { \
        $(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/testfile.nc ; \
        if [ $$? = 0 ] ; then \
           echo "PASS:     parallel run on 4 processes --------------- $$i"; \
        else \
           echo "FAILED:     parallel run on 4 processes ------------- $$i"; \
           exit 1; \
        fi ; } ; done ; } ; done

ptests: ptest4
ptest2 ptest3 ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest3 ptest4 ptest6 ptest8 ptest10

