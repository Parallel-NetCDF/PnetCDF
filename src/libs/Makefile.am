#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2815 2017-01-01 23:06:33Z wkliao $
#
# @configure_input@

lib_LTLIBRARIES = libpnetcdf.la

libpnetcdf_la_SOURCES =
libpnetcdf_la_LIBADD = ../dispatchers/libdispatchers.la ../drivers/ncmpi/libpnetcdfc.la
nodist_EXTRA_libpnetcdf_la_SOURCES = dummyc.c
if HAS_MPICXX
libpnetcdf_la_LIBADD += ../binding/cxx/libpnetcdfcxx.la
# According Automake Manual, Section 8.3.5 Libtool Convenience Libraries,
# using nodist_EXTRA to trick automake to use C++ compiler
nodist_EXTRA_libpnetcdf_la_SOURCES += dummycxx.cxx
endif
if HAS_FORTRAN
# According to Automake Manual, Section 8.14.3 Mixing Fortran 77 With C and
# C++, adding $(FLIBS) is necessary
libpnetcdf_la_LIBADD += ../binding/f77/libpnetcdff.la ../binding/f90/libpnetcdff90.la $(FLIBS) $(FCLIBS)
# nodist_EXTRA_libpnetcdf_la_SOURCES += dummyf90.f90
endif

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

