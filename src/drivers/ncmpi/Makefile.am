#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

if PNETCDF_DEBUG
   AM_CPPFLAGS = -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libpnetcdfc.la

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_HFILES = ncx_h.m4

M4_SRCS   = getput.m4 \
            i_getput.m4 \
            bput.m4 \
            ncx.m4 \
            varn.m4 \
            i_varn.m4 \
            attr.m4 \
            convert_swap.m4

H_SRCS = ncmpi_dispatch.h \
         fbits.h \
         nc.h \
         ncio.h \
         rnd.h \
         ncmpidtype.h \
         macro.h \
         utf8proc_data.h \
         utf8proc.h \
         subfile.h

C_SRCS = ncmpi_dispatch.c \
         mpinetcdf.c \
         header.c \
         mpincio.c \
         dim.c \
         error.c \
         nc.c \
         string.c \
         var.c \
         ncmpidtype.c \
         filetype.c \
         nonblocking.c   \
         malloc.c \
         utf8proc.c \
         vard.c \
         fill.c \
         util.c \
         hash_func.c

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

ncx.h: ncx_h.m4 Makefile
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

libpnetcdfc_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libpnetcdfc_la_SOURCES = ncx.h $(M4_SRCS:.m4=.c)

BUILT_SOURCES = ncx.h $(M4_SRCS:.m4=.c)

if ENABLE_SUBFILING
libpnetcdfc_la_SOURCES += subfile.c
endif

CLEANFILES = $(M4_SRCS:.m4=.c) ncx.h

EXTRA_DIST = $(M4_HFILES) \
             $(M4_SRCS) \
             subfile.c

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

