#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = pnetcdf-api
DIST_SUBDIRS = $(SUBDIRS)

EXTRA_DIST = README.CRAY \
             README.Fujitsu \
             README.INTEL \
             README.IBM \
             README.K-Computer \
             README.LINUX \
             README.SGI \
             README.SX \
             README.Ubuntu \
             README.consistency \
             README.large_files \
             README.NetCDF4.md \
             README.ADIOS.md \
             README.burst_buffering \
             pbs.script

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi
