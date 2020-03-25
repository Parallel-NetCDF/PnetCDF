#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

# adding . to the list of SUBDIRS causes the directories built before .
# This is useful for "make -j" option, as libs should be built after all
# sub-libraries
SUBDIRS = include dispatchers drivers binding packaging . libs utils
DIST_SUBDIRS = include dispatchers drivers binding packaging libs utils

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	for d in $(DIST_SUBDIRS) ; do \
	    if [ "$$d" = "." ] ; then continue ; fi ; \
	    if [ -d $$d/.deps ] ; then rmdir $$d/.deps || true ; fi ; \
	done ; \
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       if [ "$$d" = "." ] ; then \
	          continue ; \
	       fi ; \
	       rmdir $$d || true ; \
	   done ; \
	fi
