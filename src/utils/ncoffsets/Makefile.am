#
# Copyright (C) 2015, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncoffsets
EXTRA_DIST = ncoffsets.c

ncoffsets$(EXEEXT): ncoffsets.c
	$(SEQ_CC) -o $@ $<

dist_man_MANS = ncoffsets.1

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncoffsets.c
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncoffsets.c
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncoffsets.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncoffsets.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncoffsets.1

