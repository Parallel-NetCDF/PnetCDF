## Utility Programs

* [ncmpidiff](#ncmpidiff)
* [ncdiff](#ncdiff)
* [ncmpidump](#ncmpidump)
* [ncmpigen](#ncmpigen)
* [ncoffsets](#ncoffsets)
* [ncvalidator](#ncvalidator)
* [pnetcdf-config](#pnetcdf-config)
* [pnetcdf_version](#pnetcdf_version)

### ncmpidiff
* An MPI program runs in parallel to compare the contents of the two files and
  reports the first difference to the standard output.

### ncdiff
* A sequential version of `ncmpidiff`, compares the contents of the two classic
  netCDF files and reports the first difference found to the standard output.
  The classic file formats include CDF-1, CDF-2, and CDF-5.

### ncmpidump
* `ncmpidump` generates an ASCII representation of a specified netCDF file on
  standard output.

### ncmpigen
* `ncmpigen` generates either a netCDF file, or C or Fortran source code to
  create a netCDF file.

### ncoffsets
* `ncoffsets` prints the file offsets information of variables defined in a
  given netCDF file.

### ncvalidator
* `ncvalidator` checks the header of a netCDF file for whether it conforms the
  classic CDF file formats.

### pnetcdf-config
* `pnetcdf-config` displays the build and installation information of the
  PnetCDF library.

### pnetcdf_version
* `pnetcdf_version` prints the version information of PnetCDF library and the
  configure command line used to build the library

Copyright (C) 2012, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.
