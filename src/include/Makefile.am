#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2815 2017-01-01 23:06:33Z wkliao $
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

nodist_include_HEADERS = pnetcdf.h

EXTRA_DIST = pnetcdf.h.in ncconfig.h.in nctypes.h dispatch.h

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

