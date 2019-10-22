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
  + A new command-line option `-t` is added to utility program `cdfdiff` to
    compare variable differences within a tolerance. See the man page of
    `cdfdiff` for usage.

* Other updates:
  + none

* Bug fixes
  + Fix strict aliasing bug when building PnetCDF with -O3 flag.

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
  + When building PnetCDF using OpenMPI 4.0.2, some of MPI constants deprecated
    in MPI 3.0, such as MPI_COMBINER_HVECTOR_INTEGER, are still defined in
    mpi.h, even when the compatibility of MPI-1 is disable. However, they are
    defined as error messages to cause a compile-time error if user programs
    try to use them. This behavior only appears for some new versions of C
    compilers with __STDC_VERSION__ >= 201112L. When using earlier versions,
    those deprecated MPI constants are not defined at all. See discussion
    [issue 7099](https://github.com/open-mpi/ompi/issues/7099). Thanks to Carl
    Ponder who found that using gcc version 7.4.0 to build PnetCDF failed with
    error message "error: expected expression before _Static_assert" when the
    compiler sees deprecated MPI_COMBINER_HVECTOR_INTEGER. Thanks to Nick
    Papior, the solution to eliminate this error is to rebuild OpenMPI using
    the configure option "--enable-mpi1-compatibility". Note MPICH does not
    have such an issue, as all deprecated constants are still defined.

* Issues related to Darshan library:
  + none

* Clarifications
  + none

