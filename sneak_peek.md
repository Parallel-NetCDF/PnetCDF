------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + Improve the performance for the case when posting nonblocking requests of
    more than one record. See commit a00f5d2 and its comments.

* New Limitations
  + none

* Update configure options
  + Retire configure options `--enable-netcdf4` and `--enable-adios`, as there
    are already options `--with-netcdf4` and `--with-adios`. According to the
    autoconf manual, `--enable-feature` is for internal packages and
    `--with-feature` is for external.  Adding `--with-feature` is equivalent to
    enabling the feature.

* New constants
  + none

* New APIs
  + none

* API syntax changes
  + none

* API semantics updates
  + none

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + none

* New PnetCDF hint
  + none

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility program
  + none

* Other updates:
  + IBM XLF compiler on Summit at OLCF requires -qfixed compile flag when
    compiling Fortran programs written in fixed form. Thus, the fixed form flag
    detected at configure time has been added to AM_FFLAGS when compiling fixed
    form programs. Similarly, the free form flag has been added to AM_FCFLAGS
    when compiling free form programs. This issue affects only the test and
    example Fortran programs. The PnetCDF library is intact.
    See [PR #73](https://github.com/Parallel-NetCDF/PnetCDF/pull/73)
  + Add all PnetCDF I/O hints to the inquired MPI info object returned by
    ncmpi_inq_file_info(). See commit f0e65cf.

* Bug fixes
  + Fix configure bug of setting environment variable SEQ_CC to the sequential
    CC extracted from MPICC. Add configure help message for environment
    variable SEQ_CC. See commit 4978f6d. Thanks to Carl Ponder for reporting.
  + When calling MPI_File_close, fh argument fh should be a pointer.
    See 3e331a6
  + When using an MPI implementation whose Fortran feature was disabled, the
    MPI Fortran constants and datatypes may not be defined in the header file
    mpi.h. This is the case for Open MPI (tested with 4.0.2). PnetCDF used some
    Fortran datatypes without checking whether they are defined, which can fail
    at 'make' time. A fix has been added that checks whether the Fortran
    feature is disabled and wraps around the Fortran datatypes with 'ifdef
    ENABLE_FORTRAN' directive. Thanks Bert Wesarg for reporting this bug.  See
    issue [#68](https://github.com/Parallel-NetCDF/PnetCDF/issues/68) and pull
    request [#69](https://github.com/Parallel-NetCDF/PnetCDF/pull/69).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + none

* Conformity with NetCDF library
  + none

* Discrepancy from NetCDF library
  + none

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + Nonblocking APIs have yet been updated to support subfiling feature.

