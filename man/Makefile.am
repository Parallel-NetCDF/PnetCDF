#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

man_MANS = pnetcdf.3
if HAS_FORTRAN
man_MANS += pnetcdf_f77.3 pnetcdf_f90.3
endif

CLEANFILES = $(man_MANS)

EXTRA_DIST = pnetcdf.m4 pnetcdf_f90.m4

pnetcdf.3: pnetcdf.m4
	$(M4) $(M4FLAGS) -DAPI=C -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || rm -f $@

pnetcdf_f77.3: pnetcdf.m4
	$(M4) $(M4FLAGS) -DAPI=F -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || rm -f $@

pnetcdf_f90.3: pnetcdf_f90.m4
	$(M4) $(M4FLAGS) -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || rm -f $@

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/pnetcdf.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/pnetcdf.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/pnetcdf_f90.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/pnetcdf_f90.m4

