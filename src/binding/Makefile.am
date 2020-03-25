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

DIST_SUBDIRS = cxx f77 f90

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	for d in $(DIST_SUBDIRS) ; do \
	   if [ -d $$d/.deps ] ; then rmdir $$d/.deps || true ; fi ; \
	done ; \
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi
