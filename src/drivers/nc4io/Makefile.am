#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 3283 2017-07-30 21:10:11Z wkliao $
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_builddir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include
AM_CPPFLAGS += -I${top_builddir}/src/drivers/include
AM_CPPFLAGS += @NETCDF4_INC@

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libnc4io.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS = nc4io_get_put.m4

H_SRCS = nc4io_driver.h

C_SRCS = nc4io_attr.c \
         nc4io_dim.c \
         nc4io_driver.c \
         nc4io_file.c \
         nc4io_var.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libnc4io_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libnc4io_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcno

EXTRA_DIST = $(M4_HFILES) $(M4_SRCS)

tests-local: all

