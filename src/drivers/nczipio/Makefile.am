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

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libnczipio.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS =

H_SRCS = nczipio_driver.h

C_SRCS = nczipio_attr.c \
         nczipio_dim.c \
         nczipio_driver.c \
         nczipio_file.c \
         nczipio_var.c \
         nczipio_internal.c \
         nczipio_util.c \
         nczipioi_put_var.c \
         nczipioi_get_var.c \
         nczipioi_put_varn.c \
         nczipioi_iput.c \
         nczipioi_nonblocking.c \
         nczipioi_chunk.c \
         nczipioi_var.c \
         nczipioi_lists.c \
         nczipioi_nonblocking.c \
         nczipioi_wait.c \
         nczip_dummy.c

if ENABLE_ZLIB
   C_SRCS += nczip_zlib.c
endif

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libnczipio_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libnczipio_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_HFILES) $(M4_SRCS)

tests-local: all

