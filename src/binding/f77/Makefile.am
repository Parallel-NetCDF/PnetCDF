#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .a .o .c .h .f .inc

nodist_include_HEADERS = pnetcdf.inc

C2F_SRCS = issyserrf.c \
           xstrerrorf.c \
           xstrerrnof.c \
           xinq_libversf.c \
           allf.c

LIB_FSRCS = strerrorf.f \
            strerrnof.f \
            inq_libversf.f

noinst_LTLIBRARIES  = libpnetcdff.la
libpnetcdff_la_SOURCES = $(LIB_FSRCS) mpinetcdf_impl.h
nodist_libpnetcdff_la_SOURCES = $(C2F_SRCS) mpifnetcdf.h

EXTRA_DIST = buildiface \
	     defs \
	     createffiles

CLEANFILES = mpifnetcdf.h $(C2F_SRCS)

$(C2F_SRCS): mpifnetcdf.h

PNETCDF_HEADER = $(top_builddir)/src/include/pnetcdf.h

mpifnetcdf.h: defs buildiface $(PNETCDF_HEADER)
	$(RM) -f $@ allf.c
	$(srcdir)/buildiface -infile=$(PNETCDF_HEADER) -deffile=$(srcdir)/defs > allf.c
	@touch $(C2F_SRCS)

tests-local: all

dist-hook:
	-rm -rf `find $(distdir) -name .svn`
