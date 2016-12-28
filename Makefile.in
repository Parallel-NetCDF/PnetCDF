#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

srcdir = @srcdir@
VPATH = @srcdir@

include macros.make

PACKAGE_NAME      = @PACKAGE_NAME@
PACKAGE_VERSION   = @PACKAGE_VERSION@
PACKAGE_STRING    = @PACKAGE_STRING@
PACKAGE_BUGREPORT = @PACKAGE_BUGREPORT@

PACKAGE_FILENAME  = $(PACKAGE_NAME)-$(PACKAGE_VERSION)

INSTALL_SUBDIRS	= src man

SUBDIRS		= src \
		  test \
		  scripts \
		  doc \
		  man \
		  benchmarks \
		  examples

GARBAGE		= pnetcdf_pc

DIST_GARBAGE	= config.cache \
		  config.status \
		  config.log \
		  macros.make \
		  stamp-h \
                  pnetcdf.pc \
                  pnetcdf-config \
                  check_install

PACKING_LIST	= Makefile.in \
		  aclocal.m4 \
		  configure \
		  configure.in \
		  macros.make.in \
		  rules.make \
		  COPYRIGHT \
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
		  stamp-h.in \
                  pnetcdf_pc.in \
                  pnetcdf-config.in \
                  check_install.in

PACKING_SUBDIRS	= $(SUBDIRS)

all:
	$(MAKE) $(MFLAGS) -C src

TESTS_SUBDIRS = test examples benchmarks
TESTS_DIRS = $(TESTS_SUBDIRS:%=tests-%)
tests: $(TESTS_DIRS)
$(TESTS_DIRS): 
	$(MAKE) $(MFLAGS) -C $(@:tests-%=%)

check testing:
	$(MAKE) $(MFLAGS) -C test testing

verbose_check verbose_testing:
	$(MAKE) $(MFLAGS) -C test verbose_testing

PTEST_SUBDIRS = test examples benchmarks
PTEST_DIRS = $(PTEST_SUBDIRS:%=ptest-%)
ptest: tests $(PTEST_DIRS)
$(PTEST_DIRS): 
ifeq (@enable_coverage@, yes)
	echo "Parallel test is disabled because coverage analysis was enabled"
else
	$(MAKE) $(MFLAGS) -C $(@:ptest-%=%) ptest
endif
# make sure ptest runs one directory after another
ptest-examples: ptest-test
ptest-benchmarks: ptest-examples

PTESTS_DIRS = $(PTEST_SUBDIRS:%=ptests-%)
ptests: tests $(PTESTS_DIRS)
$(PTESTS_DIRS): 
ifeq (@enable_coverage@, yes)
	echo "Parallel test is disabled because coverage analysis was enabled"
else
	$(MAKE) $(MFLAGS) -C $(@:ptests-%=%) ptests
endif
# make sure ptests runs one directory after another
ptests-examples: ptests-test
ptests-benchmarks: ptests-examples

INSTALLDIRS = $(INSTALL_SUBDIRS:%=install-%)
install: $(INSTALLDIRS) install_PKGCONFIG install_CONFIG
	@echo '+----------------------------------------------------------------------------+'
	@echo '|'
	@echo '|  PnetCDF has been successfully installed under $(prefix)'
	@echo '|'
	@echo '|  * PnetCDF header files are installed in $(prefix)/include'
	@echo '|  * PnetCDF library file is installed in $(prefix)/lib'
	@echo '|  * PnetCDF utility programs are installed in $(prefix)/bin'
	@echo '|  * PnetCDF man pages are installed in $(prefix)/man/man1 and'
	@echo '|    $(prefix)/man/man3'
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

$(INSTALLDIRS): 
	$(MAKE) $(MFLAGS) -C $(@:install-%=%) install

PKGCONFIG_IN = pnetcdf_pc
PKGCONFIG    = pnetcdf.pc

install_PKGCONFIG:
	$(RM) -f $(PKGCONFIG)
	$(SED) -e "s%INSTALL_PREFIX%${prefix}%g" $(PKGCONFIG_IN) > $(PKGCONFIG)
	$(INSTALL) -d -m 755 $(LIBDIR)/pkgconfig
	$(INSTALL_DATA) $(PKGCONFIG) $(LIBDIR)/pkgconfig/$(PKGCONFIG)

install_CONFIG:
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 pnetcdf-config $(BINDIR)/pnetcdf-config

UNINSTALLDIRS = $(INSTALL_SUBDIRS:%=uninstall-%)
uninstall: $(UNINSTALLDIRS) uninstall_PKGCONFIG uninstall_CONFIG
$(UNINSTALLDIRS): 
	$(MAKE) $(MFLAGS) -C $(@:uninstall-%=%) uninstall

uninstall_PKGCONFIG:
	$(RM) -f $(LIBDIR)/pkgconfig/$(PKGCONFIG)

uninstall_CONFIG:
	$(RM) -f $(BINDIR)/pnetcdf-config

################################################################################
# Distribution:

# The following rule checks to see that were on the right system.  If we're
# not, then the distribution mightn't contain the intended ncgen(1) (not
# to mention that dvips(1) is necessary for building the C++ User's Guide).
#
check_system:	FORCE
	@case "$(OVERRIDE)" in \
	'') case `uname -sr` in \
	    'SunOS 5'*) \
		exit 0;; \
	    *)  echo 1>&2 "Error: Not on a SunOS 5 system."; \
		echo 1>&2 "Set macro OVERRIDE to override."; \
		exit 1;; \
	    esac;; \
	*) exit 0;; \
	esac

# Make a compressed, tar(1) file of the source distribution in the current 
# directory.
#
tar.Z:		check_system FORCE
	$(MAKE) $(MFLAGS) $(PACKAGE_FILENAME).tar.Z

dist: ensure_manifest MANIFEST
	id=$(PACKAGE_FILENAME) \
	&& $(RM) -fr $$id \
	&& $(LN_S) $(srcdir) $$id \
	&& tar $(TARFLAGS) - `$(SED) "s|^|$$id/|" MANIFEST` | gzip  > $(PACKAGE_FILENAME).tar.gz \
	&& tar $(TARFLAGS) - `$(SED) "s|^|$$id/|" MANIFEST` | bzip2 > $(PACKAGE_FILENAME).tar.bz2 \
	&& $(RM) -f $$id MANIFEST

MANIFEST:	FORCE
	$(MAKE) $(MFLAGS) -s MANIFEST.echo >$@

# Make a compressed, tar(1) file of the source distribution in the
# appropriate FTP directory.
#
# NOTE: Making "ftp" will cause the "tar.Z" file to be made even if it
# already exists.  This is because the MANIFEST file upon which it
# depends must be remade every time.  This is not a waste of time,
# however, if the "tar.Z" target is made in private directories and the
# "ftp" target is made in the "/upc/$(PACKAGE_FILENAME)/build/" directory.
#
ftp:		check_system FORCE
	$(MAKE) $(MFLAGS) $(FTPDIR)/$(PACKAGE_FILENAME).tar.Z \
	test -r $(FTPDIR)/$(PACKAGE_FILENAME).tar.Z || exit 0; \
	    cd $(FTPDIR) || exit 1;

$(FTPDIR)/$(PACKAGE_FILENAME).tar.Z:	$(PACKAGE_FILENAME).tar.Z
	$(RM) -f $@
	cp $(PACKAGE_FILENAME).tar.Z $@ 
	chmod u+rw,g+rw,o=r $@

# Make a compressed, tar(1) file of the binary distribution in the 
# appropriate FTP directory.
#
binftp:		FORCE
	$(MAKE) $(MFLAGS) $(FTPBINDIR)/$(PACKAGE_FILENAME).tar.Z \

ftpbin:		binftp

$(FTPBINDIR)/$(PACKAGE_FILENAME).tar.Z:
	$(RM) -f $@
	id=$(PACKAGE_FILENAME) \
	&& $(RM) -f $$id \
	&& $(LN_S) $(prefix) $$id \
	&& tar $(TARFLAGS) - README $$id/bin $$id/include \
	    $$id/lib $$id/man | compress > $@ \
	&& $(RM) -f $$id
	chmod u+rw,g+rw,o=r $@
	test -r $(FTPBINDIR)/$(PACKAGE_FILENAME).tar.Z || exit 0; \
	    cd $(FTPBINDIR) || exit 1;


# rule for generating cscope information
cscope:
	find $(srcdir) -iname "*.[ch]" > $(srcdir)/cscope.files
	( cd $(srcdir); cscope -be -i $(srcdir)/cscope.files )

# The following dependency is for configure.in and configure
# See autoconf manual 2.69, Section 4.8.5 Automatic Remaking
$(srcdir)/configure: configure.in aclocal.m4
	cd '$(srcdir)' && autoconf
     
# autoheader might not change src/lib/ncconfig.h.in so touch a stamp file.
$(srcdir)/src/lib/ncconfig.h.in: stamp-h.in
$(srcdir)/stamp-h.in: configure.in aclocal.m4
	cd '$(srcdir)' && autoheader
	echo timestamp > '$(srcdir)/stamp-h.in'
     
src/lib/ncconfig.h: stamp-h
stamp-h: src/lib/ncconfig.h.in config.status
	./config.status
     
Makefile: Makefile.in config.status
	./config.status
     
config.status: configure
	./config.status --recheck

include $(srcdir)/rules.make


