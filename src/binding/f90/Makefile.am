#
# Copyright (C) 2013, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .mod .f90 .m4 .fh

AM_DEFAULT_SOURCE_EXT = .f90

AM_FCFLAGS = -I. -I$(srcdir)

if UPPER_CASE_MOD
   PNETCDF_MOD = PNETCDF.$(FC_MODEXT)
else
   PNETCDF_MOD = pnetcdf.$(FC_MODEXT)
endif

M4FFLAGS += -DINTENTV=$(INTENTV)

M4_SRCS = getput_text.m4 \
          getput_var.m4 \
          getput_varn.m4 \
          getput_vard.m4

F90_SRCS = nf90_constants.fh \
           overloads.fh \
           visibility.fh \
           file.fh \
           dims.fh \
           attributes.fh \
           variables.fh

SRCS = $(F90_SRCS) \
       $(M4_SRCS:.m4=.fh) \
       api.fh \
       nfmpi_constants.fh

BUILT_SOURCES = $(M4_SRCS:.m4=.fh)

nodist_include_HEADERS = $(PNETCDF_MOD)
noinst_LTLIBRARIES = libf90.la
libf90_la_SOURCES = $(F90_SRCS)
nodist_libf90_la_SOURCES = pnetcdf.f90 api.fh nfmpi_constants.fh $(M4_SRCS:.m4=.fh)

.m4.fh:
	$(M4) $(AM_M4FLAGS) $(M4FFLAGS) $< >$@

pnetcdf.lo $(PNETCDF_MOD): pnetcdf.f90 $(SRCS)

EXTRA_DIST = $(M4_SRCS)

CLEANFILES = $(PNETCDF_MOD) $(M4_SRCS:.m4=.fh) *.gcda *.gcno *.gcov gmon.out

tests-local: all

