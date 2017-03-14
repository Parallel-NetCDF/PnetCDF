#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src man

if BUILD_TESTSETS
SUBDIRS += examples benchmarks test
endif

if BUILD_DOCS
SUBDIRS += doc
endif

CLEANFILES = pnetcdf_pc

DISTCLEANFILES = pnetcdf.pc

EXTRA_DIST = COPYRIGHT \
             CREDITS \
             INSTALL \
             README \
             README.CRAY \
             README.Fujitsu \
             README.INTEL \
             README.IBM \
             README.K-Computer \
             README.LINUX \
             README.SGI \
             README.SX \
             README.consistency \
             README.large_files \
             RELEASE_NOTES \
             pbs.script \
             cobalt.script \
             check_install.in \
             m4/foreach.m4 \
             m4/foreachq.m4 \
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
	( for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest $$* ; \
	done ) ;

ptests:
	( for d in $(PTEST_SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptests $$* ; \
	done ) ;

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libpnetcdf.pc

# The script shows the end users how pnetcdf is built
bin_SCRIPTS = pnetcdf-config

install-data-hook:
	@echo '+----------------------------------------------------------------------------+'
	@echo '|'
	@echo '|  PnetCDF has been successfully installed under $(prefix)'
	@echo '|'
	@echo '|  * PnetCDF header files are installed in $(prefix)/include'
	@echo '|  * PnetCDF library file is installed in $(prefix)/lib'
	@echo '|  * PnetCDF utility programs are installed in $(prefix)/bin'
	@echo '|  * PnetCDF man pages are installed in $(prefix)/share/man'
	@echo '|'
	@echo '|  To compile your PnetCDF programs, please add the following to the command'
	@echo '|  line, so the compiler can find the PnetCDF header files:'
	@echo '|      -I$(prefix)/include'
	@echo '|'
	@echo '|  Add the following line to link your program to PnetCDF library:'
	@echo '|      -L$(prefix)/lib -lpnetcdf'
	@echo '|'
	@echo '|  PnetCDF is jointly developed by a team at Northwestern University and'
	@echo '|  Argonne National Laboratory.'
	@echo '|'
	@echo '|  Visit PnetCDF web sites for more information'
	@echo '|      http://cucis.ece.northwestern.edu/projects/PnetCDF'
	@echo '|      http://trac.mcs.anl.gov/projects/parallel-netcdf'
	@echo '|'
	@echo '+----------------------------------------------------------------------------+'

.PHONY: ptest ptests tests

dist-hook:
	-rm -rf `find $(distdir) -name .svn`


