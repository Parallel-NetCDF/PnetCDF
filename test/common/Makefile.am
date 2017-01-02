#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .F90 .h

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I$(top_builddir)/src/lib

check_LTLIBRARIES = libtestutils.la

libtestutils_la_SOURCES = testutils.c testutils.h

# Must build the Fortran library separately, otherwise libtools will keep
# using FCLINK to build C library even when --disable-fortran is used
if HAS_FORTRAN
check_LTLIBRARIES += libtestutilsf.la
libtestutilsf_la_SOURCES = testutilsf.F90
libtestutils_la_LIBADD = libtestutilsf.la
CLEANFILES = testutilsf.mod
endif

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS  = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  = $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS = $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

# build check targets but not invoke
tests: $(check_LTLIBRARIES)
.PHONY: tests ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

