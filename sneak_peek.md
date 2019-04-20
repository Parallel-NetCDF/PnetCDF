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

* New/updated utility program
  + none

* Other updates:
  + none

* Bug fixes
  + When `--enable-netcdf4` is used at configure time, users may encounter
    problem during configure or make time, if the NetCDF4 library was built
    with static libraries only. Thanks Bruno Pagani for reporting. This has
    been fixed in
    [pull request #46](https://github.com/Parallel-NetCDF/PnetCDF/pull/46).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/nc4/notsupport - Test if error code NC_ENOTSUPPORT is properly
    returned when calling APIs for unsupported NetCDF-4 feature.
  + test/nc4/rec - Test creating and reading a NetCDF-4 file with 1 unlimited
    dimension. 
  + test/nc4/rec2 - Test opening a NetCDF-4 file with more than 1 unlimited
    dimensions.

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

