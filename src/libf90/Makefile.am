#
# Copyright (C) 2013, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .mod .f90 .m4

AM_DEFAULT_SOURCE_EXT = .f90

AM_FCFLAGS = -I. -I$(top_srcdir)/src/libf90

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

F90_SRCS = nf90_constants.f90 \
           overloads.f90 \
           visibility.f90 \
           file.f90 \
           dims.f90 \
           attributes.f90 \
           variables.f90

nodist_include_HEADERS = pnetcdf.mod
noinst_LTLIBRARIES = libpnetcdff90.la
libpnetcdff90_la_SOURCES = pnetcdf.f90

BUILT_SOURCES = $(M4_SRCS:.m4=.f90)

SRCS     = $(F90_SRCS) $(M4_SRCS:.m4=.f90) \
           api.f90 nfmpi_constants.f90

# libpnetcdff90_la_DEPENDENCIES = $(SRCS)
# EXTRA_libpnetcdff90_la_DEPENDENCIES = $(SRCS)

pnetcdf.o pnetcdf.mod: pnetcdf.f90 $(SRCS)
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c $(FCFLAGS_f90) pnetcdf.f90
	touch $*.mod

.m4.f90:
	$(M4) $(AM_M4FLAGS) $(M4FFLAGS) $< >$@

EXTRA_DIST = $(F90_SRCS) \
             $(M4_SRCS) \
             Makefile.in \
             pnetcdf.f90.in \
             api.f90.in \
             nfmpi_constants.f90.in

CLEANFILES = $(PNETCDF_MOD) $(M4_SRCS:.m4=.f90)
DISTCLEANFILES = pnetcdf.f90 api.f90 nfmpi_constants.f90

dist-hook:
	-rm -rf `find $(distdir) -name .svn`
