------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Accessing HDF5-based NetCDF-4 files is supported. PnetCDF can now be built
    on top of NetCDF-4, which allows PnetCDF to read and write a NetCDF-4 file.
    Note NetCDF-4 file format is currently not supported by the burst buffering
    driver.
  + Thread-safe capability is added to this release. It can be enabled with
    configure option `--enable-thread-safe`. In addition, option
    `--with-pthread` can be used to specify the path to the pthreads library.
    This feature currently only supports one-thread-per-file I/O operation.

* New optimization
  + The internal data buffering mechanism in the burst buffer driver is
    removed. This mechanism caches the request data in memory until the
    accumulated size is more than 8 MiB, so the write requests to burst buffers
    can be aligned with 8 MiB boundaries. However, experiments on Cray DataWarp
    show a negligible performance improvement unless the I/O request is small
    and fragment. On the other hand, it can degrade performance for mid- and
    large-sized requests. The burst buffer driver now writes directly to the
    burst buffers for each user write request.

* New Limitations
  + Due to a bug in HDF5 1.10.2 that fails zero-length write requests in the
    collective mode, PnetCDF is not able to support such requests when NetCDF-4
    feature is enabled.

* Update configure options
  + Enable NetCDF-4 support.
    - `--enable-netcdf4`: enable NetCDF4 format classic mode support
    - `--with-hdf5=/path/to/hdf5`: set HDF5 install path
    - `--with-netcdf4=/path/to/netcdf-4`: set NetCDF install path

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
  + doc/README.NetCDF4.md is added to describe the usage of the new feature of
    NetCDF-4 support.

* New/updated utility program
  + none

* Other updates:
  + none

* Bug fixes
  + For put and get APIs when buftype is MPI_DATATYPE_NULL, bufcount is
    ignored. This is not implemented correctly for blocking put and get APIs.
    See commit
    [403e483](https://github.com/Parallel-NetCDF/PnetCDF/commit/403e4839cdfca6175bcb177f3efa16f7d5e602d2)
  + ncmpidiff -- when comparing two files that contain record variables but
    no record has been written. See commit
    [2d2cacb](https://github.com/Parallel-NetCDF/PnetCDF/commit/2d2cacbad20b71c36a9442d9abb6c113f1838d28)
  + ncmpidiff -- when comparing two scalar variables, error NC_EBADDIM may
    mistakenly reported. See commit
    [b4d2dda](https://github.com/Parallel-NetCDF/PnetCDF/commit/b4d2dda2362e0dc0926d2723bffffea61df7006d)

* New example programs
  + examples/C/pthread.c - demonstrates the one-file-per-thread I/O example.
    When running on some parallel machines, users may need to set certain
    environment variable to enable MPI multi-threading support, for example on
    Cori @NERSC with command
    ```
    export MPICH_MAX_THREAD_SAFETY=multiple
    ```

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/error_precedence.m4 - tests the error code reporting
    precedence
  + test/nc4/put_get_all_kinds.m4 - tests all supported variable read/write
    API. Make sure they are properly wired up
  + test/nc4/interoperability_rd.m4 - tests whether NetCDF-4 file written using
    NetCDF can be read by PnetCDF
  + test/nc4/interoperability_wr.m4 - tests whether NetCDF-4 file written using
    PnetCDF can be read by NetCDF
  + test/nc4/simple_xy.c - tests reading NetCDF-4 files, borrowed the test
    program simple_xy.c from NetCDF
  + test/testcases/tst_pthread.c - tests thread-safe capability for scenario of
    each thread operating on a unique file.

* Conformity with NetCDF library
  + In contract to NetCDF-4 which allows to read/write variables in define mode
    when the file format is in NetCDF-4 format, PnetCDF still requires reading
    and writing variables in data mode.

* Conflict with NetCDF library
  + none

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + none

