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
	$(M4) $(M4FLAGS) -DAPI=C -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || $(RM) -f $@

pnetcdf_f77.3: pnetcdf.m4
	$(M4) $(M4FLAGS) -DAPI=F -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || $(RM) -f $@

pnetcdf_f90.3: pnetcdf_f90.m4
	$(M4) $(M4FLAGS) -DPNETCDF_VERSION=$(PNETCDF_VERSION) $? > $@  || $(RM) -f $@

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE2|`date '+%e %b %Y'`|g" $(distdir)/pnetcdf.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|`date '+%Y-%m-%d'`|g"  $(distdir)/pnetcdf.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE2|`date '+%e %b %Y'`|g" $(distdir)/pnetcdf_f90.m4
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|`date '+%Y-%m-%d'`|g"  $(distdir)/pnetcdf_f90.m4

