#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_srcdir)/src/drivers/include
AM_CPPFLAGS += -I$(top_builddir)/src/drivers/include
AM_CPPFLAGS += @ADIOS_INC@

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libcommon.la

H_SRCS = utf8proc_data.h \
         utf8proc.h

M4FLAGS += -I$(top_srcdir)/m4
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
            utils.c \
            error_posix2nc.c \
            hash_map.c

if ENABLE_ADIOS
   C_SRCS += error_adios2nc.c 
endif

libcommon_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libcommon_la_SOURCES = $(M4_SRCS:.m4=.c)

$(top_builddir)/src/drivers/include/ncx.h:
	set -e; cd ../include && $(MAKE) $(MFLAGS) ncx.h

# Thanks to Rafik Zurob to provide the fix and point out when using the IBM xlc
# compilers, ncx.c breaks strict ANSI C aliasing rules in regards to how
# put_ix_float / get_ix_float are used. That's why the file was miscompiled.
# The file should either be changed to follow strict ANSI C aliasing, or
# -fno-strict-aliasing should be used to compile.
if MPICC_IS_PGCC
ncx.lo: LTCOMPILE += -alias=traditional
else
if MPICC_IS_FCCPX
# Fujitsu fccpx and mpifccpx uses -Knomemalias, or -Xg -fno-strict-aliasing
ncx.lo: LTCOMPILE += -Knomemalias
else
if MPICC_IS_XLC
# IBM xlc
ncx.lo: LTCOMPILE += -qalias=noansi
else
# GNU gcc, Intel icc
ncx.lo: LTCOMPILE += -fno-strict-aliasing
endif
endif
endif

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make foo') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

EXTRA_DIST = $(M4_SRCS)

CLEANFILES = $(M4_SRCS:.m4=.c) *.gcda *.gcno *.gcov gmon.out

tests-local: all

