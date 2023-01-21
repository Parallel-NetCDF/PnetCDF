#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = pnetcdf-api
DIST_SUBDIRS = $(SUBDIRS)

EXTRA_DIST = README.CRAY.md \
             README.Fujitsu.md \
             README.INTEL.md \
             README.IBM.md \
             README.K-Computer.md \
             README.LINUX.md \
             README.SGI.md \
             README.SX.md \
             README.Ubuntu.md \
             README.consistency.md \
             README.large_files.md \
             README.NetCDF4.md \
             README.ADIOS.md \
             README.burst_buffering.md \
             pbs.script

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	@if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi
