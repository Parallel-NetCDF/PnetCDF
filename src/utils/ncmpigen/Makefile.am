#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .h

AM_DEFAULT_SOURCE_EXT = .c
AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/drivers/ncmpio

bin_PROGRAMS = ncmpigen
ncmpigen_SOURCES = main.c load.c escapes.c getfill.c init.c genlib.c ncmpigentab.c
ncmpigen_LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
ncmpigen_LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

NCMPIDUMP = ../ncmpidump/ncmpidump
dist_man_MANS = ncmpigen.1

HFILES = generic.h genlib.h ncmpigen.h ncmpigentab.h

EXTRA_DIST = $(HFILES) $(MANUAL) \
             ncmpigenyy.c Makefile.in \
             ncmpigen.l ncmpigen.y c0.cdl

CLEANFILES = c0.nc c1.cdl c1.nc c2.cdl f0.nc \
	     ctest.c ctest ctest0.nc ctest1.cdl \
	     ftest.f ftest ftest0.nc ftest1.cdl \
	     lex.ncmpi.c y.tab.c y.tab.h \
             core core.* *.gcda *.gcno *.gcov gmon.out

# generating the fortran does not work yet
# test:           $(PROGRAM) b-test c-test f-test FORCE
# test:           $(PROGRAM) b-test c-test FORCE

# Below is used if a PnetCDF developer wants to rebuild ncmpigenyy.c or
# ncmpigentab.c. In that case, configure.in at the root directory needs
# to check the availability of commands yacc/lex/bison. Otherwise the
# below is never invoked, but records how to do it.
# $(srcdir)/ncmpigentab.c \
# $(srcdir)/ncmpigentab.h:	ncmpigen.y ncmpigenyy.c ncmpigen.h
# if HAVE_YACC_LEX
# 	$(yacc) -d $(srcdir)/ncmpigen.y; \
# 	cp -f y.tab.c ncmpigentab.c; \
# 	cp -f y.tab.h ncmpigentab.h
# else
# 	@echo "Error: one of $? is modified, but cannot find bison or yacc, required to re-generate $@"
# endif

# ncmpigenyy.c:	ncmpigen.l
# if HAVE_YACC_LEX
# 	$(lex) $(srcdir)/ncmpigen.l; \
# 	cp -f lex.ncmpi.c ncmpigenyy.c
# else
# 	@echo "Error: $? is modified, but cannot find flex or lex, required to re-generate $@"
# endif

# vmstab.h \
# vmstab.c:	ncmpigen.y
# 	@echo 1>&2 "$@ is out-of-date with respect to $?"
# 	@echo 1>&2 "It must be recreated via POSIX yacc(1) on a VMS system"
# 	false
# vms_yy.c:	ncmpigenyy.c
# 	@echo 1>&2 "$@ is out-of-date with respect to $?"
# 	@echo 1>&2 "It must be recreated via POSIX lex(1) on a VMS system"
# 	false

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncmpigen.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncmpigen.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncmpigen.1

