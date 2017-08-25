#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = pnetcdf.pc

pnetcdf.pc: pnetcdf_pc
	$(SED) -e "s%INSTALL_PREFIX%${prefix}%g" $< > $@

CLEANFILES = pnetcdf.pc

.PHONY: pnetcdf.pc
