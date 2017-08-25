#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include
AM_CPPFLAGS += -I${top_builddir}/src/drivers/include

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libcommon.la

H_SRCS = utf8proc_data.h

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS   = ncx.m4 \
            convert_swap.m4

C_SRCS    = utf8proc.c \
            mem_alloc.c \
            dtype_decode.c \
            create_imaptype.c \
            error_mpi2nc.c \
            check_name.c \
            pack_unpack.c \
            utils.c

libcommon_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libcommon_la_SOURCES = $(M4_SRCS:.m4=.c)

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make foo') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

EXTRA_DIST = $(M4_SRCS)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcno

tests-local: all

