#
# Copyright (C) 2013, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(srcdir)/../common
LDADD = $(top_builddir)/src/libs/libpnetcdf.la ../common/libtestutils.la

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS = -DHAVE_DECL_MPI_OFFSET
   # AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   # AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

TESTPROGRAMS = test_subfile

check_PROGRAMS = $(TESTPROGRAMS)

TESTS = seq_runs.sh

CLEANFILES = *.nc core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = README seq_runs.sh

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

ptest ptest4: $(TESTPROGRAMS)
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i -f $(TESTOUTDIR)/testfile.nc -s 2 ; \
	) ; done

ptest2: $(TESTPROGRAMS)
	set -e ; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i -f $(TESTOUTDIR)/testfile.nc -s 2 ; \
	) ; done

ptests: ptest2 ptest4
ptest6 ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

