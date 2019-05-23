#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += @ADIOS_INC@

NCGEN = ../ncmpigen/ncmpigen

bin_PROGRAMS = ncmpidump
ncmpidump_SOURCES = ncmpidump.c vardata.c dumplib.c
ncmpidump_LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
ncmpidump_LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

dist_man_MANS = ncmpidump.1

HFILES = ncmpidump.h vardata.h dumplib.h

EXTRA_DIST = $(HFILES) $(MANUAL) test0.cdl

CLEANFILES = test0.nc test1.nc test1.cdl test2.cdl \
             core core.* *.gcda *.gcno *.gcov gmon.out

test: $(bin_PROGRAMS)
	$(NCGEN) -b -o `pwd`/test0.nc $(srcdir)/test0.cdl
	./$(bin_PROGRAMS) -n test1 test0.nc > test1.cdl
	$(NCGEN) -b -o `pwd`/test1.nc `pwd`/test1.cdl
	./$(bin_PROGRAMS) test1.nc > test2.cdl
	@cmp test1.cdl test2.cdl && \
	     echo "*** $(bin_PROGRAMS) test successful ***"

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncmpidump.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncmpidump.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncmpidump.1

