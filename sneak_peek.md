------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Suuport MPI derived datatypes that are constructed from MPI 4.0 large-count
    derived datatype constructors.
    See [PR #145](https://github.com/Parallel-NetCDF/PnetCDF/pull/145).

* New optimization
  + none

* New Limitations
  + none

* Configure options
  + The default has been changed to build both shared and static libraries.
    See [PR #143](https://github.com/Parallel-NetCDF/PnetCDF/pull/143).

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
  + none

* New example programs
  + C/flexible_bottom.c and C/vard_bottom.c - These two examples construct MPI
    derived data types using absolute memory addresses first and use MPI_BOTTOM
    when calling the PnetCDF flexible APIs.

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/flexible_large_count.c - tests flexible APIs that use MPI
    derived datatypes created by MPI large-count datatype constructors.
    See [PR #145](https://github.com/Parallel-NetCDF/PnetCDF/pull/145).

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

