#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I$(top_srcdir)/src/drivers/ncmpio
AM_CPPFLAGS += -I$(top_srcdir)/src/drivers/include
AM_CPPFLAGS += -I$(top_builddir)/src/drivers/ncmpio
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/drivers/include

bin_PROGRAMS = ncmpivalid
ncmpivalid_SOURCES = ncmpivalid.c
ncmpivalid_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpivalid.1

CLEANFILES = *.gcno

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g" $(distdir)/ncmpivalid.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|`date '+%e %b %Y'`|g"    $(distdir)/ncmpivalid.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE2|`date '+%Y-%m-%d'`|g"   $(distdir)/ncmpivalid.1

