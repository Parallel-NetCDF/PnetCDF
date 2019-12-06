------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + none

* New Limitations
  + Configure command now checks whether the supplied MPI C compiler is a
    wrapper of a C++ compiler. If this is detected, the PnetCDF configuration
    will be aborted. This check is enforced because using such an MPI C
    compiler will cause problem for linking Fortran, C and C++ programs, with
    an error message similar to this:
    ```
    conftestf.o: In function `MAIN_': conftestf.f:4: undefined reference to `sub_'
    configure:33318: error: Could not link conftestf.o and conftest.o
    ```

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
  + One Theta @ ALCF, when compiling utility programs ncoffsets and cdfdiff,
    one may encounter rrror messages below, when the compile environment is set
    to use Intel-based compilers, i.e. PrgEnv-intel
    ```
    In file included from /usr/include/inttypes.h:27:0,
    /theta-archive/intel/compilers_and_libraries_2019.5.281/linux/compiler/include/stdint.h:43:54: error: missing binary operator before token "("
       defined(__has_include_next) && __has_include_next(<stdint.h>)
                                                        ^
    ```
    This can be resolved by adding "SEQ_CC=icc" to your make command line, i.e.
    ```
    make SEQ_CC=icc
    ```

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
    in MPI 3.0, such as `MPI_COMBINER_HVECTOR_INTEGER`, are still defined in
    mpi.h, even if the compatibility option of MPI-1 is disable (default) when
    OpenMPI is configured. However, they are defined as error messages to cause
    a compile-time error if user programs try to use them. This behavior first
    appears in OpenMPI 4.0.2 when the underlying C compilers supporting 2011
    revision of the C standard (with `__STDC_VERSION__ >= 201112L`) such as GCC
    version 4.6 and later are used to build OpenMPI. When using earlier
    versions of C compilers to build OpenMPI 4.0.2, those deprecated MPI
    constants are not defined at all in mpi.h. See discussion
    [issue 7099](https://github.com/open-mpi/ompi/issues/7099). Thanks to Carl
    Ponder who found that using gcc version 7.4.0 to build PnetCDF failed with
    error message "error: expected expression before _Static_assert" when the
    compiler sees the deprecated constant `MPI_COMBINER_HVECTOR_INTEGER`.
    Thanks to Nick Papior for providing a workaround solution: to rebuild
    OpenMPI 4.0.2 and add configure option "--enable-mpi1-compatibility" when
    building PnetCDF 1.12.0 and earlier versions. Note the latest MPICH 3.3.1
    does not have such an issue, as all deprecated constants are still defined.
    This issue is now fixed in PnetCDF of release 1.12.1 which no longer
    requires the workaround build of OpenMPI.

* Issues related to Darshan library:
  + none

* Clarifications
  + The string length of I/O hint `nc_burst_buf_dirname`, the name of burst
    buffer directory must be less than the value of MPI_MAX_INFO_VAL. This is
    because all PnetCDF I/O hints were handled through MPI info mechanism and
    MPI requires the maximum string length of the value of an MPI info object
    to be MPI_MAX_INFO_VAL. If violated, error MPI_ERR_INFO_VALUE will be
    returned.

