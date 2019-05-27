#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2847 2017-01-17 07:40:16Z wkliao $
#
# @configure_input@

SUFFIXES = .o .c .h .m4

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_builddir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include
AM_CPPFLAGS += @NETCDF4_INC@ @ADIOS_INC@

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libdispatchers.la

M4FLAGS += -I${top_srcdir}/m4

M4_SRCS  = attr_getput.m4 \
           var_getput.m4

C_SRCS = lib_version.c \
         file.c \
         dimension.c \
         variable.c \
         attribute.c \
         error_codes.c

libdispatchers_la_SOURCES = $(C_SRCS)

nodist_libdispatchers_la_SOURCES = $(M4_SRCS:.m4=.c)

BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_SRCS)

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

tests-local: all

