#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS = -I$(top_builddir)/src/drivers/ncmpi -I$(top_srcdir)/src/drivers/ncmpi

bin_PROGRAMS = ncmpivalid
ncmpivalid_SOURCES = ncmpivalid.c
ncmpivalid_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

dist_man_MANS = ncmpivalid.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

