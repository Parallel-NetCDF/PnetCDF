#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src man examples benchmarks test
DIST_SUBDIRS = src man examples benchmarks test doc

if BUILD_DOCS
SUBDIRS += doc
endif

EXTRA_DIST = COPYRIGHT \
             CREDITS \
             INSTALL \
             README \
             RELEASE_NOTES \
             m4/foreach.m4 \
             m4/utils.m4

# Below is a trick to build all test executables, without running them
# tests:
# 	$(MAKE) $(AM_MAKEFLAGS) check TESTS=

# below is the alternative when AM_EXTRA_RECURSIVE_TARGETS is not defined
# AM_EXTRA_RECURSIVE_TARGETS macro was introduced into automake 1.13.
# TESTS_SUBDIRS = test examples benchmarks
# TESTS_DIRS = $(TESTS_SUBDIRS:%=tests-%)
# tests: all $(TESTS_DIRS)
# $(TESTS_DIRS):
# 	$(MAKE) $(MFLAGS) -C $(@:tests-%=%) tests

if BUILD_COVERAGE
   PTEST_SUBDIRS =
else
   PTEST_SUBDIRS = test examples benchmarks
endif

ptest:
	for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest $$* || exit 1 ; \
	done

ptests:
	for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptests $$* || exit 1 ; \
	done

tests-local:
	@echo '|'
	@echo '|  All test programs have been successfully built.'
	@echo '|'

check-local:
	@echo '|'
	@echo '|  All sequential test programs have run successfully.'
	@echo '|'

install-data-hook:
	@echo '+----------------------------------------------------------------------------+'
	@echo '|'
	@echo '|  PnetCDF has been successfully installed under'
	@echo '|          $(DESTDIR)$(prefix)'
	@echo '|'
	@echo '|  * PnetCDF header files have been installed in'
	@echo '|          $(DESTDIR)$(prefix)/include'
	@echo '|  * PnetCDF library files have been installed in'
	@echo '|          $(DESTDIR)$(exec_prefix)/lib'
	@echo '|  * PnetCDF utility programs have been installed in'
	@echo '|          $(DESTDIR)$(exec_prefix)/bin'
	@echo '|  * PnetCDF man pages have been installed in'
	@echo '|          $(DESTDIR)$(prefix)/share/man'
	@echo '|'
	@echo '|  To compile your PnetCDF programs, please add the following to the command'
	@echo '|  line, so the compiler can find the PnetCDF header files:'
	@echo '|      -I$(DESTDIR)$(prefix)/include'
	@echo '|'
	@if test "x$(has_fortran)" = xyes -a 'x$(FC_MODINC)' != 'x-I' ; then \
	 echo '|  Add the following line to compile your Fortran programs' ; \
	 echo '|      $(FC_MODINC)$(DESTDIR)$(prefix)/include' ; \
	 echo '|' ; \
	fi
	@echo '|  Add the following line to link your program to PnetCDF library:'
	@echo '|      -L$(DESTDIR)$(exec_prefix)/lib -lpnetcdf'
	@echo '|'
	@if test "x$(enable_shared)" = xyes ; then \
	 echo '|  Add the following to your run-time environment variable LD_LIBRARY_PATH,' ; \
	 echo '|  when linking your executable with the PnetCDF shared libraries.' ; \
         if test "x$(enable_netcdf4)" = xyes -a "x$(enable_adios)" = xyes ; then \
         echo '|      $(DESTDIR)$(exec_prefix)/lib:$(netcdf4_libdir):$(adios_libdir)' ; \
         elif test "x$(enable_netcdf4)" = xyes ; then \
         echo '|      $(DESTDIR)$(exec_prefix)/lib:$(netcdf4_libdir)' ; \
         elif test "x$(enable_adios)" = xyes ; then \
         echo '|      $(DESTDIR)$(exec_prefix)/lib:$(adios_libdir)' ; \
         else \
	 echo '|      $(DESTDIR)$(exec_prefix)/lib' ; \
	 fi ; \
	 echo '|' ; \
	fi
	@echo '|'
	@echo '|  PnetCDF is jointly developed by a team at Northwestern University and'
	@echo '|  Argonne National Laboratory.'
	@echo '|'
	@echo '|  Visit PnetCDF project web site for more information'
	@echo '|      https://parallel-netcdf.github.io'
	@echo '|'
	@echo '+----------------------------------------------------------------------------+'

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi

.PHONY: ptest ptests tests

dist-hook:
	$(SED_I) -e "1,10s|_PNETCDF_RELEASE_DATE_|@PNETCDF_RELEASE_DATE@|" $(distdir)/RELEASE_NOTES
	$(SED_I) -e "1,10s|_PNETCDF_VERSION_|$(PNETCDF_VERSION)|" $(distdir)/RELEASE_NOTES
	$(SED_I) -e "s|_PNETCDF_VERSION_|$(PNETCDF_VERSION)|" $(distdir)/INSTALL
	-rm -rf `find $(distdir) -type d -name .svn`


