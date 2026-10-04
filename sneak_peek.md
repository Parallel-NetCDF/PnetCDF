------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New feature
  + none

* New optimization
  + none

* New I/O driver
  + none

* New Limitations
  + none

* New configure options
  + none

* Configure option updates:
  + none

* New constants
  + none

* New APIs
  + none

* APIs deprecated
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
  + Fix PnetCDF package file, pnetcdf.pc.in, by removing GIO library as a
    required library. Because GIO's source codes are included in all PnetCDF
    official releases, it is not necessary to make gio as a required package.
    Thanks Victor Eijkhout for reporting in
    [Issue #245](https://github.com/Parallel-NetCDF/PnetCDF/issues/245). See
    the fix in [PR #246](https://github.com/Parallel-NetCDF/PnetCDF/pull/246).
  + Fix a compilation error when using Intel oneAPI compilers causing some
    Fortran APIs not visible in the shared library. Thanks to Xylar Asay-Davis
    for reporting and suggesting the fix in
    [Issue #242](https://github.com/Parallel-NetCDF/PnetCDF/issues/242). See
    the fix in [PR #244](https://github.com/Parallel-NetCDF/PnetCDF/pull/244).
  + Fix a compilation error when using C23 compilers, where boolean true and
    false are keywords rather than macros. Thanks to Xylar Asay-Davis for
    reporting and suggesting the fix in
    [Issue #240](https://github.com/Parallel-NetCDF/PnetCDF/issues/240). See
    the fix in [PR #241](https://github.com/Parallel-NetCDF/PnetCDF/pull/241).

* New example programs
  + none

* New I/O benchmarks
  + none

* New test programs
  + none

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

* Clarifications about of PnetCDF hints
  + none

