#
# Copyright (C) 2019, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_builddir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include
AM_CPPFLAGS += -I${top_builddir}/src/drivers/include
AM_CPPFLAGS += -I${srcdir}/adios_headers
AM_CPPFLAGS += -I${srcdir}/adios_headers/core
AM_CPPFLAGS += -I${srcdir}/adios_headers/core/transforms
AM_CPPFLAGS += -I${srcdir}/adios_headers/transforms
AM_CPPFLAGS += @ADIOS_INC@

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libncadios.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_SRCS = ncadios_convert.m4

H_SRCS = ncadios_driver.h \
         ncadios_internal.h

C_SRCS = ncadios_attr.c \
         ncadios_dim.c \
         ncadios_driver.c \
         ncadios_file.c \
         ncadios_var.c \
         ncadios_misc.c \
         ncadios_internal.c \
         ncadios_lists.c \
         ncadios_bp2ncd.c \
         ncadios_sync.c \
         ncadios_nonblocking.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libncadios_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libncadios_la_SOURCES = $(M4_SRCS:.m4=.c)

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make target') in a clean tree if it depends on a built source."
BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out

EXTRA_DIST = $(M4_SRCS) adios_headers

tests-local: all

