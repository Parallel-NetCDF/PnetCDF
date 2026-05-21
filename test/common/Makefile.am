#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .F90 .h

AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utils/ncmpidiff

LDADD = @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@ -lm

check_LTLIBRARIES = libtestutils.la

libtestutils_la_LIBADD = $(top_builddir)/src/utils/ncmpidiff/libncmpidiff_core.la

libtestutils_la_SOURCES = testutils.c testutils.h

CLEANFILES = core.* *.gcda *.gcno *.gcov gmon.out

# Must build the Fortran library separately, otherwise libtools will keep
# using FCLINK to build C library even when --disable-fortran is used
if HAS_FORTRAN
check_LTLIBRARIES += libtestutilsf.la
libtestutilsf_la_SOURCES = testutilsf.F90
libtestutils_la_LIBADD += libtestutilsf.la

CLEANFILES += testutilsf.mod
endif

AM_FFLAGS  =
AM_FCFLAGS =

if RELAX_COORD_BOUND
   AM_FFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
   AM_FCFLAGS += $(FC_DEFINE)RELAX_COORD_BOUND
endif

if DECL_GET_ENVIRONMENT_VARIABLE
   AM_FCFLAGS += $(FC_DEFINE)HAS_GET_ENVIRONMENT_VARIABLE
endif

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will be used by both C and Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   # AM_CFLAGS  = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

# overwrite PPFCCOMPILE and LTPPFCCOMPILE to avoid using DEFS
# IBM xlf compiler uses -WF,-D not -D
# FCDEFS=`echo $(DEFS) | $(SED) 's/-D/-WF,-D/g'`
# PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
#         $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
# LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
#         $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
#         $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
#         $(AM_FCFLAGS) $(FCFLAGS)

$(top_builddir)/src/utils/ncmpidiff/libncmpidiff_core.la:
	set -e; cd $(top_builddir)/src/utils/ncmpidiff && $(MAKE) $(MFLAGS)

# build check targets but not invoke
tests-local: all $(check_LTLIBRARIES)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

