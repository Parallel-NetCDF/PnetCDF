#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

lib_LTLIBRARIES  = libpnetcdf.la
nodist_include_HEADERS = pnetcdf.h

M4FLAGS += -I${top_srcdir}/m4
if ENABLE_ERANGE_FILL
M4FLAGS += -DERANGE_FILL
endif

M4_HFILES = ncx_h.m4

M4_SRCS   = getput.m4 \
            i_getput.m4 \
            bput.m4 \
            ncx.m4 \
            m_getput_varx.m4 \
            varn.m4 \
            i_varn.m4 \
            attr.m4 \
            convert_swap.m4

H_SRCS = fbits.h \
         nc.h \
         nctypes.h \
         ncio.h \
         rnd.h \
         ncmpidtype.h \
         macro.h \
         utf8proc_data.h \
         utf8proc.h \
         subfile.h

C_SRCS = mpinetcdf.c \
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

libpnetcdf_la_SOURCES = $(C_SRCS) $(H_SRCS)
nodist_libpnetcdf_la_SOURCES = ncx.h $(M4_SRCS:.m4=.c)

BUILT_SOURCES = ncx.h $(M4_SRCS:.m4=.c)

if ENABLE_SUBFILING
libpnetcdf_la_SOURCES += subfile.c
endif

libpnetcdf_la_LIBADD =
if HAS_MPICXX
libpnetcdf_la_LIBADD += ../libcxx/libpnetcdfcxx.la
endif
if HAS_FORTRAN
libpnetcdf_la_LIBADD += ../libf/libpnetcdff.la ../libf90/libpnetcdff90.la $(FLIBS) $(FCLIBS)
endif

CLEANFILES = $(M4_SRCS:.m4=.c) ncx.h
# DISTCLEANFILES = ncconfig.h pnetcdf.h

EXTRA_DIST = $(M4_HFILES) \
             $(M4_SRCS) \
             pnetcdf.h.in \
             ncconfig.h.in \
             subfile.c

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

