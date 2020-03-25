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

if ENABLE_BURST_BUFFER
   SUBDIRS += ncbbio
endif

if ENABLE_ADIOS
   SUBDIRS += ncadios
endif

DIST_SUBDIRS = include common ncmpio ncfoo ncbbio nc4io ncadios

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
