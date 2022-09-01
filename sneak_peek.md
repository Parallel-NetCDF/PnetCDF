------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Flexible APIs now can be used as high-level APIs, when argument bufcount
    is -1 and buftype is an MPI predefined data type. For example,
    ```
    ncmpi_put_vara_all(ncid, varid, start, count, buf, -1, MPI_FLOAT);
    ```
    is equivalent to
    ```
    ncmpi_put_vara_float_all(ncid, varid, start, count, buf);
    ```
    See [PR #82](https://github.com/Parallel-NetCDF/PnetCDF/pull/82).

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

* Updated utility program
  + ncvalidator - Add printing of the dimension size of a variable when its
    size is larger than the limitation allowed by the file format. See commit
    5584d44.

* Other updates:
  + In all prior versions, the file name was checked whether it contains
    character ':'. The prefix name ending with ':' is considered by ROMIO as
    the file system type name. The prefix name, if found, is then stripped, so
    the file name can be used in the successive POSIX function calls. However,
    the prefix was not checked against the file system type names recognized
    by ROMIO. Starting from this release, the prefix is checked against the
    known file system type names to ROMIO. If the prefix is not one of the
    recognized types, e.g.  "ufs", "nfs", "xfs", "pvfs2", "gpfs", "panfs",
    "lustre", "daos", "testfs", "ime", or "quobyte", then the prefix name is
    not stripped. This change is for in case when the file name contains ':',
    but it is not for specifying the file system type.
    See [PR #79](https://github.com/Parallel-NetCDF/PnetCDF/pull/79).

* Bug fixes
  + none

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/test_get_varn.c - test get_varn API. See PR #90.
  + test/testcases/flexible_var.c - test flexible var API
  + test/testcases/flexible_api.f - test flexible API when bufcount == -1

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

