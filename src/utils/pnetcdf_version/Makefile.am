#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = pnetcdf_version
EXTRA_DIST = pnetcdf_version.c

XDEFS = -DMPICC="\"$(MPICC)\""   -DCFLAGS="\"$(CFLAGS)\""
XDEFS += -DCONFIGURE_ARGS_CLEAN="\"$(CONFIGURE_ARGS_CLEAN)\""
XDEFS += -DPNETCDF_VERSION="\"$(PNETCDF_VERSION)\""
XDEFS += -DPNETCDF_RELEASE_DATE="\"$(PNETCDF_RELEASE_DATE)\""

if HAS_MPICXX
XDEFS += -DMPICXX="\"$(MPICXX)\"" -DCXXFLAGS="\"$(CXXFLAGS)\""
endif

if HAS_FORTRAN
XDEFS += -DMPIF77="\"$(MPIF77)\"" -DF77FLAGS="\"$(F77FLAGS)\""
XDEFS += -DMPIF90="\"$(MPIF90)\"" -DF90FLAGS="\"$(F90FLAGS)\""
endif

pnetcdf_version$(EXEEXT): pnetcdf_version.c
	$(SEQ_CC) $(XDEFS) -o $@ $<

dist_man_MANS = pnetcdf_version.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

