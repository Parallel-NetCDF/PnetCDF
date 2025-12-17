------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + none

* New Limitations
  + none

* Configure options

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
  + Fix data movement when new record variables are added to an existing file
    that does not change the starting offset of record variable section.
    See [PR #199](https://github.com/Parallel-NetCDF/PnetCDF/pull/199).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/tst_grow_data.c -- adding new variables by re-entering the
    define mode multiple time, but does not cause file header extent to grow.
    It also tests a case when adding a new record variable that does not change
    the starting offset of the record variable section in the file.

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

