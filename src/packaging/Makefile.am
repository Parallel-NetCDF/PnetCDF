#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pnetcdf.pc

# GNU autoconf allows DESTDIR variable when running 'make install' which will
# prepend it before all installation names.
install-data-hook:
	$(SED_I) -e 's|INSTALL_PREFIX|$(prefix)|g ; s|INSTALL_EXEC_PREFIX|$(exec_prefix)|g' $(DESTDIR)$(libdir)/pkgconfig/pnetcdf.pc


