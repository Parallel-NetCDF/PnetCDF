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
  + Fix strict aliasing bug when building PnetCDF with -O3 flag. See
    [a40aa5f](https://github.com/Parallel-NetCDF/PnetCDF/commit/a40aa5f73938ba1298f92ad471b3e3578ef8dbda)

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
    mpi.h, even if the compatibility option of MPI-1 is disable (default) when
    OpenMPI is configured. However, they are defined as error messages to cause
    a compile-time error if user programs try to use them. This behavior first
    appears in OpenMPI 4.0.2 when newer versions of C compilers with
    __STDC_VERSION__ >= 201112L is used to build OpenMPI. When using earlier
    versions of C compilers, those deprecated MPI constants are not defined at
    all. See discussion
    [issue 7099](https://github.com/open-mpi/ompi/issues/7099). Thanks to Carl
    Ponder who found that using gcc version 7.4.0 to build PnetCDF failed with
    error message "error: expected expression before _Static_assert" when the
    compiler sees the deprecated constant MPI_COMBINER_HVECTOR_INTEGER. Thanks
    to Nick Papior of providing a workaround solution: to rebuild OpenMPI 4.0.2
    and add configure option "--enable-mpi1-compatibility" when building
    PnetCDF 1.12.0 and erarlier versions. Note the latest MPICH 3.3.1 does not
    have such an issue, as all deprecated constants are still defined. This
    issue is now fixed in PnetCDF of release 1.12.1 which no longer requires
    the workaround build of OpenMPI.

* Issues related to Darshan library:
  + none

* Clarifications
  + none

