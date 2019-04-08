------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + When inserting nonblocking requests into the PnetCDF internal pending
    queues, the queues are now kept sorted (using an insert sort) into an
    increasing order of variable starting file offsets. This can avoid the
    quick sort when flushing the pending requests. See [pull request #37]
    (https://github.com/Parallel-NetCDF/PnetCDF/pull/37). To avoid internal
    sorts completely, users are recommended to post nonblocking requests in
    the increasing order of variable IDs and fixed-size variables first
    followed by record variables.

* New Limitations
  + When building with NetCDF-4 feature enabled, using a NetCDF-4 library that
    has already been built with PnetCDF enabled, i.e. --enable-pnetcdf, is not
    supported. See [Issue #33](https://github.com/Parallel-NetCDF/PnetCDF/issues/33).

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
  + none

* Other updates:
  + Add a check whether the MPI library is built with shared-library support.
    If not and `--enable-shared` is used, the configure process of PnetCDF will
    fail.
  + In the NetCDF-4 driver, `nc4io_inq_var()` adds a NULL-argument check for
    `no_fill` and `fill_value`. If both arguments are NULL, it skips the call
    to `nc_inq_var_fill`.
  + File header extent area between end of header and first variable will be
    padded with null bytes if PnetCDF is configured with option
    `--enable-null-byte-header-padding`.

* Bug fixes
  + Fix ncmpidiff when comparing dimension names of 2 variables between files
    whose dimension define orders are different. See
    [Issue #42](https://github.com/Parallel-NetCDF/PnetCDF/pull/42).
  + Fix error checking for programs in examples/C to ignore NC_ENOTENABLED
    if PnetCDF was not built with --enable-profiling. Thanks to Bruno Pagani
    and see [Issue #34](https://github.com/Parallel-NetCDF/PnetCDF/issues/34).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/burst_buffer/varn.c -- to test varn API when burst buffer driver is
    used. The test includes cases when argument counts is NULL or some elements
    in counts are NULL.

* Conformity with NetCDF library
  + none

* Discrepancy from NetCDF library
  + none

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + Padding -- NetCDF classic file format specification states "Header padding
    uses null (\x00) bytes. In data, padding uses variable's fill value."
    PnetCDF implements the header padding specification but only enforces it
    when the configure option `--enable-null-byte-header-padding` is set. Note
    PnetCDF has not yet implemented the padding for data section.


