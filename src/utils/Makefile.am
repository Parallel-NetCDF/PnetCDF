#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = ncmpigen ncmpidump ncmpidiff ncmpivalid pnetcdf_version ncoffsets
DIST_SUBDIRS = $(SUBDIRS)

# The script shows the end users how pnetcdf is built
bin_SCRIPTS = pnetcdf-config

EXTRA_DIST = pnetcdf-config.in

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/" $(distdir)/pnetcdf-config.in

