------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + When file header extent size grows, moving the data section to a higher
    file offset has changed to be done in chunks of 16 MB per process.
    See [PR #174](https://github.com/Parallel-NetCDF/PnetCDF/pull/174),

* New Limitations
  + none

* Configure options
  + For PnetCDF developers, the requirement for libtool version has been
    changed to 2.5.4, due to an issue on Mac OS when using OpenMPI. See
    [Issue #155](https://github.com/Parallel-NetCDF/PnetCDF/issues/155),
    [Issue #163](https://github.com/Parallel-NetCDF/PnetCDF/issues/163),
    and [PR #164](https://github.com/Parallel-NetCDF/PnetCDF/pull/164).

* Configure updates:
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

* New PnetCDF hints
  + none

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + none

* Other updates:
  + none

* Bug fixes
  + Fix setting of user hint nc_ibuf_size.
    See [PR #161](https://github.com/Parallel-NetCDF/PnetCDF/pull/161).

* New example programs
  + none

* New programs for I/O benchmarks
  + WRF-IO contains an extraction of the I/O kernel of WRF (Wether Research
    and Forecast Model, a weather prediction computer simulation program
    developed at NCAR) that can be used to evaluate the file write performance
    of WRF. It's data partitioning pattern is a 2D block-block checkerboard
    pattern, along the longitude and latitude.
    See [PR #165](https://github.com/Parallel-NetCDF/PnetCDF/pull/165).

* New test program
  + test/testcases/tst_grow_header.c tests header extent growth by re-entering
    the define mode multiple times and add more fix-sized and record variables.

* Issues with NetCDF library
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

