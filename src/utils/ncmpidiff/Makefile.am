#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

bin_PROGRAMS = ncmpidiff
ncmpidiff_SOURCES = ncmpidiff.c
ncmpidiff_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpidiff.1

CLEANFILES = *.gcno

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncmpidiff.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/ncmpidiff.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/ncmpidiff.1

tests-local: all

