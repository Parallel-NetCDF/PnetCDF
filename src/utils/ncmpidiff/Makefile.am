#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncmpidiff
ncmpidiff_SOURCES = ncmpidiff.c
ncmpidiff_CPPFLAGS = -I$(top_builddir)/src/lib
ncmpidiff_LDADD = $(top_builddir)/src/lib/libpnetcdf.la

dist_man_MANS = ncmpidiff.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

