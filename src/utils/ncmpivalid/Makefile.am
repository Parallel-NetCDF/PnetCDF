#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS = -I$(top_srcdir)/src/lib

bin_PROGRAMS = ncmpivalid
ncmpivalid_SOURCES = ncmpivalid.c
ncmpivalid_LDADD = $(top_builddir)/src/lib/libpnetcdf.la

dist_man_MANS = ncmpivalid.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

