------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + NetCDF-4 driver -- Accessing HDF5-based NetCDF-4 files is now supported.
    PnetCDF can be built on top of NetCDF-4 library to let users to use PnetCDF
    APIs to read and write a NetCDF-4 file. Users now can add NC_NETCDF4 flag
    when calling ncmpi_create() to create NetCDF-4 files. For opening NetCDF-4
    files, no additional flag is needed, as PnetCDF automatically detects the
    file format and uses the HDF5 I/O driver underneath. This feature is
    provided for convenience purpose. The parallel I/O performance to NetCDF-4
    files is expected no difference from using NetCDF-4 library directly.
  + Per-file thread-safe capability is added. This feature can be enabled at
    configure time by adding command-line option `--enable-thread-safe`. In
    addition, option `--with-pthread` can be used to specify the install path
    to the pthreads library. This feature currently only supports
    one-thread-per-file I/O operations and the classic CDF-1, 2, and 5 files.
  + ADIOS driver -- Read ADIOS 1.x BP formated file. 
    ADIOS_READ_METHOD_BP must be set when open BP file.
    Does not support low-level and non blocking API.

* New optimization
  + On some systems, e.g. Cori @NERSC, collective MPI-IO may perform poorly
    when the I/O buffer is noncontiguous, compared to a contiguous one. To
    avoid this, `ncmpi_wait()` and `ncmpi_wait_all()` now check whether the
    buffer is noncontiguous and size is less than 16 MiB. If both are true, a
    temporary contiguous buffer is allocated to copy the data over and used in
    the MPI read or write calls. See
    [PR #26](https://github.com/Parallel-NetCDF/PnetCDF/pull/26). Programs
    developed to test this issue is available in
    https://github.com/Parallel-NetCDF/E3SM-IO/tree/master/mpi_io_test

* New Limitations
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

* Update configure options
  + Enable NetCDF-4 support.
    - `--enable-netcdf4`: enable NetCDF4 format classic mode support
    - `--with-netcdf4=/path/to/netcdf-4`: path to NetCDF-4 library installation
  + Enable ADIOS support.
    - `--enable-adios`: enable NetCDF4 format classic mode support
    - `--with-adios=/path/to/netcdf-4`: path to NetCDF-4 library installation
  + Enable multi-threading support.
    - `--enable-thread-safe`: enable per-file thread-safe support
    - `--with-pthread`: path to the pthread library installation

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
  + The automatic file layout alignment for fixed-size variables is disabled.
    This is because modern MPI-IO implementations have already aligned the file
    access with the file lock boundaries and the automatic alignment can create
    a file view with "holes" in between variables, which can adversely degrade
    I/O performance.
  + The internal data buffering mechanism used in the burst buffer driver is
    removed. This mechanism caches the request data in memory until the
    accumulated size is more than 8 MiB, so the write requests to burst buffers
    can be aligned with 8 MiB boundaries. However, experiments on Cray DataWarp
    show a negligible performance improvement unless the I/O request is small
    and fragment. On the other hand, it can degrade performance for mid- and
    large-sized requests. The burst buffer driver now writes directly to the
    burst buffers for each user write request.

* Bug fixes
  + Fix bug of checking interleaved requests for scalar variables. See
    [PR #27](https://github.com/Parallel-NetCDF/PnetCDF/pull/27).
  + When building PnetCDF using the IBM xlc compiler with -O optimization
    option on Little Endian platforms, users may encounter errors related to
    strict ANSI C aliasing rules. Thanks to Jim Edwards for reporting and Rafik
    Zurob for providing the fix. See
    [Issue #23](https://github.com/Parallel-NetCDF/PnetCDF/issues/23) and
    [Pull Request #24](https://github.com/Parallel-NetCDF/PnetCDF/issues/24).
  + Shell ksh has a different way to redirect stdout and stderr from bash.
    PnetCDF configure.ac and acinclude.m4 have been developed mainly on bash.
    This bug can cause configure command to fail when using ksh. Thanks to
    @poohRui for reporting the bug. See
    [Issue #21](https://github.com/Parallel-NetCDF/PnetCDF/issues/21) and
    [PR #22](https://github.com/Parallel-NetCDF/PnetCDF/pull/22).
    However, running configure under ksh is still buggy. A GNU automake bug
    report of hanging problem can be found in
    https://lists.gnu.org/archive/html/bug-automake/2015-04/msg00000.html
    PnetCDF users are recommended to run configure under other shells.
  + For put and get APIs when buftype is MPI_DATATYPE_NULL, bufcount is
    ignored. This is not implemented correctly in blocking put and get APIs.
    See bug fix committed on Aug. 25, 2018.
  + ncmpidiff -- when comparing two files that contain record variables but
    no record has been written. See bug fix committed on Aug. 25, 2018.
  + ncmpidiff -- when comparing two scalar variables, error NC_EBADDIM may
    mistakenly reported. See bug fix committed on Aug. 12, 2018.
  + When the MPI communicator used in ncmpi_create or ncmpi_open is freed by
    the user after the call and before file is closed, programs would crash at
    ncmpi_close with MPI error of "Invalid communicator". The fix moves the
    duplication of MPI communicator to the place before calling driver create
    and open subroutines. See bug fix committed on Jul 21, 2018.

* New example programs
  + examples/C/time_var.c and examples/F77/time_var.f - show how to define,
    write, and read record variables.
  + examples/C/pthread.c - demonstrates the one-file-per-thread I/O example.
    When running on some parallel machines, users may need to set certain
    environment variable to enable MPI multi-threading support, for example on
    Cori @NERSC with command
    ```
    export MPICH_MAX_THREAD_SAFETY=multiple
    ```
  + examples/C/transpose2D.c - a 2D version of examples/C/transpose.c
  + examples/adios/read_all.c - Dump all metadata in a ADIOS BP file.

* New programs for I/O benchmarks
  + none

* New test program
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
  + test/adios/vars.c - tests if PnetCDF access ADIOS variables with stride.
  + test/adios/atts.c - tests if PnetCDF access ADIOS attributes.

* Conformity with NetCDF library
  + none

* Discrepancy from NetCDF library
  + In contract to NetCDF-4 which allows to read/write variables in define mode
    when the file format is in NetCDF-4 format, PnetCDF still requires reading
    and writing variables in data mode.
  + In contrast to the semantics of nc_set_fill() defined in NetCDF-4,
    ncmpi_set_fill() changes the fill mode of all variables newly defined in
    the current scope of defined mode. Variables affected include the ones
    defined before and after the call to ncmpi_set_fill(). Note this API has no
    effect on the already existing variables created in the previous define
    mode. This behavior follows the convention adopted by NetCDF-3. To change
    fill mode for individual variables after the call to ncmpi_set_fill(), API
    ncmpi_def_var_fill() can be used for this purpose. Refer NetCDF 4.1.3 user
    guide for semantics of
    [nc_set_fill()](https://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf-c/nc_005fset_005ffill.html).
    A discussion with NetCDF developers regarding this issue can be found in
    [1114](https://github.com/Unidata/netcdf-c/pull/1114).
  + The error code return precedence can be different between NetCDF and
    PnetCDF in some cases. A test program for error code return precedence is
    available in test/testcases/error_precedence.m4. This program can be used
    to test both PnetCDF and NetCDF libraries. Note when testing NetCDF
    programs, because NetCDF does not follow the same precedence, failures are
    expected. A discussion with NetCDF developers regarding this issue can be
    found in [334](https://github.com/Unidata/netcdf-c/issues/334).

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + PnetCDF currently does not support Fortran default integer type set to 8
    bytes (for GNU Fortran compiler, this change of default setting is done by
    using compile option -fdefault-integer-8). Checking this has been added
    and configure command will fail, once default 8-byte integer is detected.

