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

RELEASE_DATE = $(PNETCDF_RELEASE_DATE2)
RELEASE_STR  = "$(PNETCDF_VERSION) of $(PNETCDF_RELEASE_DATE)"

pnetcdf.3: pnetcdf.m4
	$(M4) $(M4FLAGS) -DAPI=C -DRELEASE_STR=$(RELEASE_STR) -DRELEASE_DATE=$(RELEASE_DATE) $? > $@  || $(RM) -f $@

pnetcdf_f77.3: pnetcdf.m4
	$(M4) $(M4FLAGS) -DAPI=F -DRELEASE_STR=$(RELEASE_STR) -DRELEASE_DATE=$(RELEASE_DATE) $? > $@  || $(RM) -f $@

pnetcdf_f90.3: pnetcdf_f90.m4
	$(M4) $(M4FLAGS) -DRELEASE_STR=$(RELEASE_STR) -DRELEASE_DATE=$(RELEASE_DATE) $? > $@  || $(RM) -f $@

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

