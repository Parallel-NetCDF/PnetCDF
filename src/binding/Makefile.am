#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2749 2016-12-28 17:51:37Z wkliao $
#
# @configure_input@

SUBDIRS =
if HAS_MPICXX
SUBDIRS += cxx
endif
if HAS_FORTRAN
SUBDIRS += f77 f90
endif

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

