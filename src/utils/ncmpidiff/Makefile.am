#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

bin_PROGRAMS = ncmpidiff cdfdiff

noinst_LTLIBRARIES = libncmpidiff_core.la
libncmpidiff_core_la_SOURCES = ncmpidiff_core.c

ncmpidiff_SOURCES = ncmpidiff.c
ncmpidiff_LDADD  = $(noinst_LTLIBRARIES) $(top_builddir)/src/libs/libpnetcdf.la
ncmpidiff_LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

cdfdiff_SOURCES = cdfdiff.c

cdfdiff$(EXEEXT): cdfdiff.c Makefile
	$(SEQ_CC) $(SEQ_CFLAGS) -I$(top_srcdir)/src/utils/ncvalidator -o $@ $< $(SEQ_LDFLAGS) $(SEQ_LIBS)

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpidiff.1 cdfdiff.1

EXTRA_DIST = ncmpidiff_core.h

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncmpidiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncmpidiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncmpidiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/cdfdiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/cdfdiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/cdfdiff.1
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/cdfdiff.c
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/cdfdiff.c


tests-local: all

