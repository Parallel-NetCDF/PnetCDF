#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS =
if HAS_MPICXX
SUBDIRS += libcxx
endif
if HAS_FORTRAN
SUBDIRS += libf libf90
endif

SUBDIRS += lib
SUBDIRS += utils

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

