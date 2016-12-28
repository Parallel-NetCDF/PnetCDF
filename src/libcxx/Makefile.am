#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .cpp .hpp

CXX_HEADER = pnetcdf

nodist_include_HEADERS = $(CXX_HEADER)

noinst_LTLIBRARIES = libpnetcdfcxx.la
# libpnetcdfcxx_la_CPPFLAGS = -I$(top_srcdir)/src/libcxx

# Note the file order is important
HPP_FILES = ncmpi_notyet.h \
            ncmpiCheck.h \
            ncmpiType.h \
            ncmpiAtt.h \
            ncmpiEnumType.h \
            ncmpiGroupAtt.h \
            ncmpiGroup.h \
            ncmpiByte.h \
            ncmpiChar.h \
            ncmpiCompoundType.h \
            ncmpiDouble.h \
            ncmpiFloat.h \
            ncmpiInt.h \
            ncmpiInt64.h \
            ncmpiOpaqueType.h \
            ncmpiShort.h \
            ncmpiUbyte.h \
            ncmpiUint.h \
            ncmpiUint64.h \
            ncmpiUshort.h \
            ncmpiException.h \
            ncmpiDim.h \
            ncmpiFile.h \
            ncmpiVarAtt.h \
            ncmpiVar.h \
            ncmpiVlenType.h

LIB_SRCS   = $(HPP_FILES:.h=.cpp)
libpnetcdfcxx_la_SOURCES = $(LIB_SRCS) $(HPP_FILES)
nodist_libpnetcdfcxx_la_SOURCES = $(CXX_HEADER)
BUILT_SOURCES = $(CXX_HEADER)

EXTRA_DIST = pnetcdf.in

CLEANFILES = $(CXX_HEADER)

$(CXX_HEADER): $(HPP_FILES) Makefile
	cat $(srcdir)/pnetcdf.in >$@
	set -e; for i in $(HPP_FILES); do ( \
	sed -e '/#ifndef/d' -e '/#define/d' -e '/#endif/d' -e '/#include/d' $(srcdir)/$$i >> $@ ; \
	) ; done

dist-hook:
	-rm -rf `find $(distdir) -name .svn`
