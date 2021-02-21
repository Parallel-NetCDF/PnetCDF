------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + none

* New Limitations
  + none

* Update configure options
  + none

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
  + none

* Bug fixes
  + When using an MPI compiler whose Fortran feature was disabled, the MPI
    Fortran constants and datatypes may not be defined in the header file
    mpi.h. This is the case for OpenMPI (tested with 4.0.2). PnetCDF used some
    Fortran datatypes without checking whether they are defined, which can fail
    at 'make' time. A fix has been added that checks whether the Fortran
    feature is disabled and wraps around the Fortran datatypes with 'ifdef
    ENABLE_FORTRAN' directive. Thanks Bert Wesarg for reporting this bug.
    See issue [#68](https://github.com/Parallel-NetCDF/PnetCDF/issues/68) and
    pull request [#69](https://github.com/Parallel-NetCDF/PnetCDF/pull/69).

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
  + none

