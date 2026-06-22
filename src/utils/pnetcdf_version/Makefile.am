#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

AM_CFLAGS = @C_VISIBILITY_FLAGS@

bin_PROGRAMS = pnetcdf_version

XDEFS = -DMPICC="\"$(MPICC)\""   -DCFLAGS="\"$(CFLAGS)\""
XDEFS += -DCONFIGURE_ARGS_CLEAN="\"$(CONFIGURE_ARGS_CLEAN)\""
XDEFS += -DPNETCDF_VERSION="\"$(PNETCDF_VERSION)\""
XDEFS += -DMPI_VERSION="\"@MPI_VERSION@\""
XDEFS += -DMPI_VENDOR="\"@MPI_VENDOR@\""
XDEFS += -DMPI_VENDOR_VERSION="\"@MPI_VENDOR_VERSION@\""

XDEFS += -DMPICC_BASE="\"@MPICC_BASE@\"" -DMPICC_BASE_VERSION="\"@MPICC_BASE_VERSION@\""

if HAS_MPICXX
XDEFS += -DMPICXX="\"$(MPICXX)\"" -DCXXFLAGS="\"$(CXXFLAGS)\""
XDEFS += -DMPICXX_BASE="\"@MPICXX_BASE@\"" -DMPICXX_BASE_VERSION="\"@MPICXX_BASE_VERSION@\""
endif

if HAS_FORTRAN
XDEFS += -DMPIF77="\"$(MPIF77)\"" -DFFLAGS="\"$(FFLAGS)\""
XDEFS += -DMPIF90="\"$(MPIF90)\"" -DFCFLAGS="\"$(FCFLAGS)\""
XDEFS += -DMPIF77_BASE="\"@MPIF77_BASE@\"" -DMPIF77_BASE_VERSION="\"@MPIF77_BASE_VERSION@\""
XDEFS += -DMPIF90_BASE="\"@MPIF90_BASE@\"" -DMPIF90_BASE_VERSION="\"@MPIF90_BASE_VERSION@\""
endif

pnetcdf_version$(EXEEXT): pnetcdf_version.c
	$(SEQ_CC) $(SEQ_CFLAGS) $(XDEFS) -o $@ $< $(SEQ_LDFLAGS) $(SEQ_LIBS)

dist_man_MANS = pnetcdf_version.1

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/pnetcdf_version.c
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/pnetcdf_version.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/pnetcdf_version.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/pnetcdf_version.1

