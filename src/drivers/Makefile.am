#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2749 2016-12-28 17:51:37Z wkliao $
#
# @configure_input@

SUBDIRS = include common ncmpio

if BUILD_DRIVER_FOO
   SUBDIRS += ncfoo
endif

if ENABLE_NETCDF4
   SUBDIRS += nc4io
endif

if BUILD_DRIVER_BB
   SUBDIRS += ncbbio
endif

DIST_SUBDIRS = include common ncmpio ncfoo ncbbio nc4io

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi
