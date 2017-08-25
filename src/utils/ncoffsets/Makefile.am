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
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncoffsets.c
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/ncoffsets.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/ncoffsets.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/ncoffsets.1

