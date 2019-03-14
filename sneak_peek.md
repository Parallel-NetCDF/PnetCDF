------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
<<<<<<< HEAD
<<<<<<< HEAD
  + NetCDF-4 driver -- Accessing HDF5-based NetCDF-4 files is now supported.
    PnetCDF can be built on top of NetCDF-4 library to allow PnetCDF to read
    and write a NetCDF-4 file. Users now can add NC_NETCDF4 flag to create
    NetCDF-4 files. For opening NetCDF-4 files, no additional flag is needed,
    as PnetCDF can automatically detect the file format.
  + Per-file thread-safe capability is added. This feature can be enabled by
    adding command-line option `--enable-thread-safe` at configure time. In
    addition, option `--with-pthread` can be used to specify the path to the
    pthreads library. This feature currently only supports one-thread-per-file
    I/O operations.
  + ADIOS driver -- Read ADIOS 1.x BP formated file. 
    ADIOS_READ_METHOD_BP must be set when open BP file.
    Does not support low-level and non blocking API.
=======
  + ADIOS driver -- Reading ADIOS 1.x BP formated file. 
    PnetCDF will detected BP file format automatically. No additional flag needed when opening a BP formated file.
    CUrrently, the ADIOS driver does not support low-level API.
>>>>>>> cf04fe9... sneakpeak
=======
  + ADIOS BP read capability -- Reading ADIOS 1.x BP formated file. 
    PnetCDF will detected BP file format automatically. No additional flag needed when opening a BP formated file.
    Currently, the ADIOS driver does not support low-level API.
>>>>>>> 97dfbef... sneakpeak

* New optimization
  + none

* New Limitations
<<<<<<< HEAD
<<<<<<< HEAD
  + For creating new files, the NetCDF-4 driver in PnetCDF supports only the
    classic model I/O operations. Advanced NetCDF-4 features, such as chunking,
    compression, etc. are not supported. This is due to the unavailability of
    the APIs for those operations.
  + The burst buffering driver does not support NetCDF-4 file formats.
  + Due to a bug in HDF5 1.10.2 that fails zero-length write requests to record
    variables in the collective mode, PnetCDF is not able to support such
    requests when NetCDF-4 feature is enabled. New HDF5 releases are expected
    to contain the fix. See discussion in https://github.com/NCAR/ParallelIO/pull/1304
    The bug fix will appear in HDF5 1.10.4 release.
  + ADIOS driver are ready only. No vard, varn, low-level, and nonblocking support.
=======
  + ADIOS driver is ready only. There is also no vard, varn and low-level API support.
>>>>>>> cf04fe9... sneakpeak
=======
  + ADIOS driver is ready only. For more detail on limitations, please refer to /doc/README.ADIOS.md
>>>>>>> 97dfbef... sneakpeak

* Update configure options
  + Enable NetCDF-4 support.
    - `--enable-netcdf4`: enable NetCDF4 format classic mode support
    - `--with-netcdf4=/path/to/netcdf-4`: path to NetCDF-4 library installation
  + Enable ADIOS support.
<<<<<<< HEAD
    - `--enable-adios`: enable NetCDF4 format classic mode support
    - `--with-adios=/path/to/netcdf-4`: path to NetCDF-4 library installation
  + Enable multi-threading support.
    - `--enable-thread-safe`: enable per-file thread-safe support
    - `--with-pthread`: path to the pthread library installation
=======
    - `--enable-adios`: enable reading BP formated file
    - `--with-adios=/path/to/netcdf-4`: path to ADIOS library installation
  + When building with NetCDF-4 feature, using NetCDF-4 library built with
    PnetCDF enabled, i.e. --enable-pnetcdf, is not supported. See
    [Issue #33](https://github.com/Parallel-NetCDF/PnetCDF/issues/33).
>>>>>>> cf04fe9... sneakpeak

* New constants
  + NC_BP and nf_bp is the flag indicating BP file access mode.
  + NC_FORMAT_BP and nf_format_bp indicates BP file format.

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
  + NC_EADIOS, NF_EADIOS
    ADIOS library internal error that does not corresponds to any PnetCDF error code

* New PnetCDF hint
  + none

* New run-time environment variables
  + none

* Build recipes
  + none

* New/updated utility program
  + ncmpidump
    It can now be used to dump BP formated file.

* Other updates:
  + none

* Bug fixes
  + When `--enable-netcdf4` is used at configure time, users may encounter
    problem during configure or make time, if the NetCDF4 library was built
    with static libraries only. Thanks Bruno Pagani for reporting. This has
    been fixed in
    [pull request #46](https://github.com/Parallel-NetCDF/PnetCDF/pull/46).

* New example programs
<<<<<<< HEAD
<<<<<<< HEAD
  + examples/C/pthread.c - demonstrates the one-file-per-thread I/O example.
    When running on some parallel machines, users may need to set certain
    environment variable to enable MPI multi-threading support, for example on
    Cori @NERSC with command
    ```
    export MPICH_MAX_THREAD_SAFETY=multiple
    ```
  + examples/C/transpose2D.c - a 2D version of examples/C/transpose.c
  + examples/adios/read_all.c - Dump all metadata in a ADIOS BP file.
=======
  + examples/adios/read_meta.c - Dump all metadata in a ADIOS BP file.
=======
  + examples/adios/read_metadata.c - Dump all metadata in a ADIOS BP file.
>>>>>>> 73b9061... add exmaple put_var_nb
  + examples/adios/read_var.c -  Read a variable form a BP file generated by the 
    ADIOS example program examples/C/arrays/arrays_write.
  + examples/adios/read_var_nb.c -  Read a variable form a BP file generated by the 
    ADIOS example program examples/C/arrays/arrays_write using non-blocking API.
>>>>>>> cf04fe9... sneakpeak

* New programs for I/O benchmarks
  + none

* New test program
  + test/nc4/notsupport - Test if error code NC_ENOTSUPPORT is properly
    returned when calling APIs for unsupported NetCDF-4 feature.
  + test/nc4/rec - Test creating and reading a NetCDF-4 file with 1 unlimited
    dimension. 
  + test/nc4/rec2 - Test opening a NetCDF-4 file with more than 1 unlimited
    dimensions.
  + test/F90/test_fill.f90 - another test for bug fix r3730.
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
  + test/testcases/tst_free_comm.c - free MPI communicator right after calling
    ncmpi_create to see if PnetCDF duplicates the communicator correctly.
  + test/adios/open.c - tests if PnetCDF recognize ADIOS file.
  + test/adios/header.c - tests if PnetCDF can parse ADIOS header.
  + test/adios/var.c - tests if PnetCDF can access ADIOS variables.
  + test/adios/varm.c - tests if PnetCDF can access ADIOS variables with discontiguous memory.
<<<<<<< HEAD
  + test/adios/vars.c - tests if PnetCDF access ADIOS variables with stride.
  + test/adios/atts.c - tests if PnetCDF access ADIOS attributes.
=======
  + test/adios/vars.c - tests if PnetCDF can access ADIOS variables with stride.
  + test/adios/ivar.c - tests if PnetCDF can access ADIOS variables using non-blocking API.
  + test/adios/ivarm.c - tests if PnetCDF can access ADIOS variables with discontiguous memory using non-blocking API.
  + test/adios/ivars.c - tests if PnetCDF can access ADIOS variables with stride using non-blocking API.
  + test/adios/att.c - tests if PnetCDF can access ADIOS attributes.
  + test/adios/indep.c - tests if PnetCDF can access ADIOS variables independently.
  + test/burst_buffer/varn.c -- to test varn API when burst buffer driver is
    used. The test includes cases when argument counts is NULL or some of the
    elements in counts are NULL.
>>>>>>> cf04fe9... sneakpeak

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

