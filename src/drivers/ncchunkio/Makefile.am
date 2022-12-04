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

if PNETCDF_PROFILING
   AM_CPPFLAGS += -DPNETCDF_PROFILING
endif

noinst_LTLIBRARIES = libncchkio.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS = ncchkioi_profile.m4 \
          ncchkioi_convert.m4

M4H_SRCS = ncchkioi_profile.m4h

H_SRCS = ncchkio_driver.h 

C_SRCS = ncchkio_attr.c \
         ncchkio_dim.c \
         ncchkio_driver.c \
         ncchkio_file.c \
         ncchkio_var.c \
         ncchkio_internal.c \
         ncchkioi_util.c \
         ncchkioi_put_var.c \
         ncchkioi_get_var.c \
         ncchkioi_put_varn.c \
         ncchkioi_get_varn.c \
         ncchkioi_iput_cb.c \
         ncchkioi_iget_cb.c \
         ncchkioi_iput.c \
         ncchkioi_iget.c \
         ncchkioi_nonblocking.c \
         ncchkioi_cache.c \
         ncchkioi_chunk.c \
         ncchkioi_chunk_size.c \
         ncchkioi_chunk_owner.c \
         ncchkioi_var_init.c \
         ncchkioi_var_resize.c \
         ncchkioi_var_wr.c \
         ncchkioi_var_rd.c \
         ncchkioi_lists.c \
         ncchkioi_nonblocking.c \
         ncchkioi_wait.c \
         ncchk_filter_dummy.c 

if ENABLE_ZLIB
   C_SRCS += ncchk_filter_zlib.c
endif

if ENABLE_SZ
   C_SRCS += ncchk_filter_sz.c
endif

$(M4_SRCS:.m4=.c): Makefile
$(M4H_SRCS:.m4h=.h): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

.m4h.h:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libncchkio_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libncchkio_la_SOURCES = $(M4_SRCS:.m4=.c) $(M4H_SRCS:.m4h=.h)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c) $(M4H_SRCS:.m4h=.h)

CLEANFILES = $(M4_SRCS:.m4=.c) $(M4H_SRCS:.m4h=.h) core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_HFILES) $(M4_SRCS) $(M4H_SRCS)

tests-local: all

