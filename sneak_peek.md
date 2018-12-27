------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + When inserting nonblocking requests into pending queues, keep the queues
    sorted (insert sort) in the increasing order of variable IDs. This can
    avoid a sorting (quick sort) when flushing the pending requests.

* New Limitations
  + When building with NetCDF-4 feature, using NetCDF-4 library built with
    PnetCDF enabled, i.e. --enable-pnetcdf, is not supported. See
    [Issue #33](https://github.com/Parallel-NetCDF/PnetCDF/issues/33).

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

* New/updated utility program
  + none

* Other updates:
  + none

* Bug fixes
  + Fix error checking for programs in examples/C to ignore NC_ENOTENABLED
    if PnetCDF was not built with --enable-profiling. Thanks to Bruno Pagani
    and see [Issue #34](https://github.com/Parallel-NetCDF/PnetCDF/issues/34).

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

