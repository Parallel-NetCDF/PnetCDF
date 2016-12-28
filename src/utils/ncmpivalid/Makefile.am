#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncmpivalid
ncmpivalid_SOURCES = ncmpivalid.c
ncmpivalid_CPPFLAGS = -I$(top_builddir)/src/lib -I$(top_srcdir)/src/lib
ncmpivalid_LDADD = $(top_builddir)/src/lib/libpnetcdf.la

dist_man_MANS = ncmpivalid.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

