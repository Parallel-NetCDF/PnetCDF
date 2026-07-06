#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2815 2017-01-01 23:06:33Z wkliao $
#
# @configure_input@

include_HEADERS = pnetcdf_commux.h

nodist_include_HEADERS = pnetcdf.h pnetcdf_comm.h

EXTRA_DIST = nctypes.h dispatch.h pnc_debug.h pnetcdf_comm.h.in

dist-hook:
	$(SED_I) -e "s|DIST_DATE|@PNETCDF_RELEASE_DATE@|g" $(distdir)/pnetcdf.h.in
