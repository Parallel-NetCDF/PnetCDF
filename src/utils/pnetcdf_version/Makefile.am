#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = pnetcdf_version

XDEFS = -DMPICC="\"$(MPICC)\""   -DCFLAGS="\"$(CFLAGS)\""
XDEFS += -DCONFIGURE_ARGS_CLEAN="\"$(CONFIGURE_ARGS_CLEAN)\""
XDEFS += -DPNETCDF_VERSION="\"$(PNETCDF_VERSION)\""

if HAS_MPICXX
XDEFS += -DMPICXX="\"$(MPICXX)\"" -DCXXFLAGS="\"$(CXXFLAGS)\""
endif

if HAS_FORTRAN
XDEFS += -DMPIF77="\"$(MPIF77)\"" -DFFLAGS="\"$(FFLAGS)\""
XDEFS += -DMPIF90="\"$(MPIF90)\"" -DFCFLAGS="\"$(FCFLAGS)\""
endif

pnetcdf_version$(EXEEXT): pnetcdf_version.c
	$(SEQ_CC) $(XDEFS) -o $@ $<

dist_man_MANS = pnetcdf_version.1

tests-local: all

dist-hook:
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/pnetcdf_version.c
	$(SED) -i "s/PNETCDF_RELEASE_VERSION/$(PNETCDF_VERSION)/" $(distdir)/pnetcdf_version.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE/`date '+%e %b %Y'`/"    $(distdir)/pnetcdf_version.1
	$(SED) -i "s/PNETCDF_RELEASE_DATE2/`date '+%Y-%m-%d'`/"   $(distdir)/pnetcdf_version.1

