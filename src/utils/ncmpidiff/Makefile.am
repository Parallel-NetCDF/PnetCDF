#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncmpidiff
ncmpidiff_SOURCES = ncmpidiff.c
ncmpidiff_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

dist_man_MANS = ncmpidiff.1

CLEANFILES = *.gcno

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

tests-local: all

