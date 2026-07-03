#
# Copyright (C) 2015, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CFLAGS = @C_VISIBILITY_FLAGS@

bin_PROGRAMS = ncoffsets
EXTRA_DIST = ncoffsets.c

ncoffsets$(EXEEXT): ncoffsets.c
	$(SEQ_CC) $(SEQ_CFLAGS) -o $@ $< $(SEQ_LDFLAGS) $(SEQ_LIBS)

dist_man_MANS = ncoffsets.1

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncoffsets.c
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncoffsets.c
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncoffsets.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncoffsets.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncoffsets.1
if IS_GIT_CLONE
	last_commit_date=`git -C $(abs_top_srcdir) log -1 --format='%cd' --date=short $(abs_srcdir)/ncoffsets.1` ; \
	$(SED_I) -e "s|GIT_SRC_LAST_COMMIT|$$last_commit_date|g"                $(distdir)/ncoffsets.1
endif


