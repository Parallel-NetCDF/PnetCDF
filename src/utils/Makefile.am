#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUBDIRS = ncmpigen ncmpidump ncmpidiff ncvalidator pnetcdf_version ncoffsets
DIST_SUBDIRS = ncmpigen ncmpidump ncmpidiff ncvalidator pnetcdf_version ncoffsets ncmpilogdump

if ENABLE_BURST_BUFFER
SUBDIRS += ncmpilogdump
endif

# The script shows the end users how pnetcdf is built
bin_SCRIPTS = pnetcdf-config

EXTRA_DIST = pnetcdf-config.in

# For VPATH build (parallel build), try delete all sub-directories
distclean-local:
	for d in $(DIST_SUBDIRS) ; do \
	   if [ -d $$d/.deps ] ; then rmdir $$d/.deps || true ; fi ; \
	done ; \
	if [ "$(abs_builddir)" != "$(abs_srcdir)" ] ; then \
	   for d in $(DIST_SUBDIRS) ; do \
	       rmdir $$d || true ; \
	   done ; \
	fi

# Note sed option -i for in-place substitution is not portable
dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g" $(distdir)/pnetcdf-config.in

install-exec-hook:
	$(SED_I) -e 's|INSTALL_PREFIX|$(prefix)|g ; s|INSTALL_EXEC_PREFIX|$(exec_prefix)|g' $(DESTDIR)$(bindir)/pnetcdf-config
	chmod +x $(DESTDIR)$(bindir)/pnetcdf-config

