------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Flexible APIs now can be used as high-level APIs, when argument bufcount
    is NC_COUNT_IGNORE and buftype is an MPI predefined data type. See
    [PR #82](https://github.com/Parallel-NetCDF/PnetCDF/pull/82). Below is an
    example of writing from a memory buffer of type float.
    ```
    ncmpi_put_vara_all(ncid, varid, start, count, buf, NC_COUNT_IGNORE, MPI_FLOAT);
    ```
    is equivalent to
    ```
    ncmpi_put_vara_float_all(ncid, varid, start, count, buf);
    ```

* New optimization
  + none

* New Limitations
  + none

* Update configure options
  + `--enable-install-examples` to install example programs under folder
    `${prefix}/pnetcdf_examples` along with run script files. An example is
    `${prefix}/pnetcdf_examples/C/run_c_examples.sh`. The default of this
    option is `disabled`.

* New constants
  + NC_COUNT_IGNORE - This is used in flexible APIs. When argument bufcount is
    NC_COUNT_IGNORE, buftype must be a predefine MPI datatype and the APIs
    operate as the high-level APIs. Fortran equivalents are NF_COUNT_IGNORE and
    NF90_COUNT_IGNORE.

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
  + ncvalidator - When the file size is larger and smaller than expected, the
    file may still be a valid netCDF file.
    * When larger, this can happen if opening an existing file that contains no
      variable. Deleting a global attribute reduces the file header size. The
      file is still a valid netCDF file. PR #99 detects this mismatch and
      truncates the file size to the header size.
    * When smaller, this can happen if the last variable is partially written.
      The expected file size is calculated based on the full sizes of all
      variables. The file is still a valid netCDF file.
    * In the former case, PR #99 changes ncvalidator to report a warning,
      rather than an error.
  + ncvalidator - Add printing of the dimension size of a variable when its
    size is larger than the limitation allowed by the file format. See commit
    5584d44.
  + Add file src/utils/README.md which gives short descriptions of the utility
    programs and collapsible bullets to display their manual pages.

* Other updates:
  + Upgrade autotool version requirement to autoconf 2.71, automake 1.16.5, and
    libtool 2.4.6.
    See [PR #95](https://github.com/Parallel-NetCDF/PnetCDF/pull/95)
    Thanks to Blaise Bourdin for pointing out in
    [Issue #94](https://github.com/Parallel-NetCDF/PnetCDF/issues/94)
    that configure failed when using Intel OneAPI 2022.2.0 compilers. The fix
    is to use autoconf 2.70 and newer.
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
  + Fix ncmpi_inq_num_rec_vars and ncmpi_inq_num_fix_vars when opening an
    existing file. See PR #103.
  + ncmpidiff -  when checking the dimensions defined in the second files
    whether they are defined in the first file. See 88cd9c1.

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/tst_symlink.c - test NC_CLOBBER on a symbolic link.
  + test/testcases/tst_del_attr.c - test delete attributes. See PR #99.
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

