------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + NetCDF-4 driver -- Accessing HDF5-based NetCDF-4 files is now supported.
    PnetCDF can be built on top of NetCDF-4, which allows PnetCDF to read and
    write a NetCDF-4 file. However, the burst buffering currently do not
    support NetCDF-4 file formats.
  + Thread-safe capability is added to this release. It can be enabled by
    command-line option `--enable-thread-safe` at configure time. In addition,
    option `--with-pthread` can be used to specify the path to the pthreads
    library. This feature currently only supports one-thread-per-file I/O
    operation.

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
  + Due to a bug in HDF5 1.10.2 that fails zero-length write requests to record
    variables in the collective mode, PnetCDF is not able to support such
    requests when NetCDF-4 feature is enabled. New HDF5 releases are expected
    to contain the fix.

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
  + Shell ksh has a different way to redirect stdout and stderr from bash.
    PnetCDF configure.ac and acinclude.m4 have been developed mainly on bash.
    This bug can cause configure command to fail when using ksh. Thanks to
    @poohRui for reporting the bug. See Issue #21 and PR #22.
  + For put and get APIs when buftype is MPI_DATATYPE_NULL, bufcount is
    ignored. This is not implemented correctly for blocking put and get APIs.
    See bug fix committed on Aug. 25, 2018.
  + ncmpidiff -- when comparing two files that contain record variables but
    no record has been written. See bug fix committed on Aug. 25, 2018.
  + ncmpidiff -- when comparing two scalar variables, error NC_EBADDIM may
    mistakenly reported. See bug fix committed on Aug. 12, 2018.

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
  + test/nc4/tst_zero_req.c - tests a HDF5 1.10.2 bug that causes test program
    to hang when writing to and reading back a 2D record variable in collective
    mode with some of the processes making zero-length requests.
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

* Discrepancy from NetCDF library
  + In contrast to nc_set_fill() in NetCDF, ncmpi_set_fill() changes the fill
    mode of all variables newly defined in the current scope of defined mode.
    Variables affected include the ones created before and after the call to
    ncmpi_set_fill(). Note this API has no effect on the existing variables
    defined in the previous define mode. This behavior follows the convention
    adopted by NetCDF-3, but not NetCDF-4. To change fill mode for individual
    variables after the call to ncmpi_set_fill(), API ncmpi_def_var_fill() can
    be used for this purpose. Reference of NetCDF 4.1.3 user guide for
    (nc_set_fill())[https://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf-c/nc_005fset_005ffill.html]
  + The error code precedence can be different between NetCDF and PnetCDF.

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + none

