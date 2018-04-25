#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

bin_PROGRAMS = ncmpilogdump
nodist_ncmpilogdump_SOURCES = ncmpilogdump.c
ncmpilogdump_LDADD = $(top_builddir)/src/libs/libpnetcdf.la

M4_SRCS  = ncmpilogdump.m4

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

BUILT_SOURCES = $(M4_SRCS:.m4=.c)

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

EXTRA_DIST = $(M4_SRCS)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out

tests-local: all

