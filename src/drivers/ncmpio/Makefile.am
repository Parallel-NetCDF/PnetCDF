#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
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

noinst_LTLIBRARIES = libncmpio.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS   = ncmpio_getput.m4 \
            ncmpio_i_getput.m4 \
            ncmpio_varn.m4 \
            ncmpio_i_varn.m4 \
            ncmpio_attr.m4

H_SRCS = ncmpio_driver.h \
         ncmpio_NC.h \
         ncmpio_subfile.h

C_SRCS = ncmpio_driver.c \
         ncmpio_create.c \
         ncmpio_open.c \
         ncmpio_enddef.c \
         ncmpio_close.c \
         ncmpio_sync.c \
         ncmpio_file_misc.c \
         ncmpio_header_get.c \
         ncmpio_header_put.c \
         ncmpio_dim.c \
         ncmpio_var.c \
         ncmpio_bput.c \
         ncmpio_filetype.c \
         ncmpio_wait.c   \
         ncmpio_vard.c \
         ncmpio_fill.c \
         ncmpio_util.c \
         ncmpio_hash_func.c \
         ncmpio_file_io.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libncmpio_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libncmpio_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make foo') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

if ENABLE_SUBFILING
libncmpio_la_SOURCES += ncmpio_subfile.c
endif

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_HFILES) \
             $(M4_SRCS) \
             ncmpio_subfile.c

tests-local: all

