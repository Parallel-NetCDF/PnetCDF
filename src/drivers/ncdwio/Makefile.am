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

noinst_LTLIBRARIES = libncdwio.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS = #ncdwio_log_flush.m4

H_SRCS = ncdwio_driver.h

C_SRCS = ncdwio_attr.c \
         ncdwio_dim.c \
         ncdwio_driver.c \
         ncdwio_file.c \
         ncdwio_var.c \
		 ncdwio_log.c \
		 ncdwio_mem.c \
         ncdwio_nonblocking.c \
		 ncdwio_util.c \
		 ncdwio_log_flush.c \
		 ncdwio_log_put.c \
		 ncdwio_sharedfile.c \
		 ncdwio_bufferedfile.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libncdwio_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libncdwio_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcno

EXTRA_DIST = $(M4_HFILES) $(M4_SRCS)

tests-local: all

