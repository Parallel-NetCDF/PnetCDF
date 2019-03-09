#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .cpp .hpp

AM_CPPFLAGS = -I$(top_builddir)/src/include

if SEEK_SET_REDEFINED
AM_CPPFLAGS += -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
endif

if REPLACE_CXX_FUNC_MACRO
AM_CPPFLAGS += -D__func__=__FUNCTION__
endif

CXX_HEADER = pnetcdf

nodist_include_HEADERS = $(CXX_HEADER)

noinst_LTLIBRARIES = libcxx.la

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
libcxx_la_SOURCES = $(LIB_SRCS) $(HPP_FILES)
nodist_libcxx_la_SOURCES = $(CXX_HEADER)
BUILT_SOURCES = $(CXX_HEADER)

EXTRA_DIST = pnetcdf.in

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out \
             $(CXX_HEADER) a.out

$(CXX_HEADER): $(HPP_FILES) Makefile
	cat $(srcdir)/pnetcdf.in >$@
	set -e; for i in $(HPP_FILES); do ( \
	$(SED) -e '/#ifndef/d' -e '/#define/d' -e '/#endif/d' -e '/#include/d' $(srcdir)/$$i >> $@ ; \
	) ; done

tests-local: all

