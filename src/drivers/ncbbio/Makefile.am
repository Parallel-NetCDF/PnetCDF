#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
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

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libncbbio.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS = #ncbbio_log_flush.m4

H_SRCS = ncbbio_driver.h

C_SRCS = ncbbio_attr.c \
         ncbbio_dim.c \
         ncbbio_driver.c \
         ncbbio_file.c \
         ncbbio_var.c \
		 ncbbio_log.c \
		 ncbbio_mem.c \
         ncbbio_nonblocking.c \
		 ncbbio_util.c \
		 ncbbio_log_flush.c \
		 ncbbio_log_put.c \
		 ncbbio_sharedfile.c \
		 ncbbio_misc.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libncbbio_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libncbbio_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcno

EXTRA_DIST = $(M4_HFILES) $(M4_SRCS)

tests-local: all

