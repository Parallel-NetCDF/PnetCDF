#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2847 2017-01-17 07:40:16Z wkliao $
#
# @configure_input@

SUFFIXES = .o .c .h .m4

noinst_LTLIBRARIES = libdispatchers.la

M4FLAGS += -I${top_srcdir}/m4

M4_SRCS  = attribute.m4

C_SRCS = file.c \
         dimension.c \
         error_codes.c \
         variable.c

libdispatchers_la_SOURCES = $(C_SRCS)

nodist_libdispatchers_la_SOURCES = $(M4_SRCS:.m4=.c)

BUILT_SOURCES = $(M4_SRCS:.m4=.c)

CLEANFILES = $(M4_SRCS:.m4=.c)

EXTRA_DIST = $(M4_SRCS)

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

