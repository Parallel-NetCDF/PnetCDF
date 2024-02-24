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

* Configure options
  + `--disable-file-sync` is now deprecated. This configure option alone does
    not provide a sufficient data consistency. Users are suggested to call
    `ncmpi_sync` and `MPI_Barrier` to achieve a desired consistency.
  + `--enable-install-examples` to install example programs under folder
    `${prefix}/pnetcdf_examples` along with run script files. An example is
    `${prefix}/pnetcdf_examples/C/run_c_examples.sh`. The default of this
    option is `disabled`.
  + Add three new environment variables `SEQ_CFLAGS`, `SEQ_LDFLAGS` and
    `SEQ_LIBS` for setting the compile, link, and library flags, respectively
    to be used to build the sequential utility programs, i.e. `cdfdiff`,
    `ncoffsets`, `ncvalidator`, and `pnetcdf_version`.
    See [PR #122](https://github.com/Parallel-NetCDF/PnetCDF/pull/122)

* Configure updates:
  + Upgrade config.guess config.sub to 2024-01-01.
    See [PR #116](https://github.com/Parallel-NetCDF/PnetCDF/pull/116)
  + FLASH-IO benchmark - add compile flag "-fallow-argument-mismatch" for GNU
    Fortran 10 and later.
    See [PR #114](https://github.com/Parallel-NetCDF/PnetCDF/pull/114)
  + Handle the case when MPICC env is not set and "--with-mpi" is not used.
    See commit 6142135.
  + Upgrade autotools version requirement to autoconf 2.71, automake 1.16.5, and
    libtool 2.4.6.
    See [PR #95](https://github.com/Parallel-NetCDF/PnetCDF/pull/95)
    Thanks to Blaise Bourdin for pointing out in
    [Issue #94](https://github.com/Parallel-NetCDF/PnetCDF/issues/94)
    that configure failed when using Intel OneAPI 2022.2.0 compilers. The fix
    is to use autoconf 2.70 and newer.

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
  + File open flag NC_SHARE is now deprecated. It is still defined, but takes
    no effect.
  + NC_SHARE alone is not sufficient to provide data consistency for accessing
    a shared file in parallel and thus is now deprecated.  Because PnetCDF
    follows the MPI file consistency, which only addresses the case when all
    file accesses are relative to a specific file handle created from a
    collective open, NC_SHARE becomes invalid. See doc/README.consistency.md
    for more information.

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + none

* New PnetCDF hint
  + nc_header_collective -- to instruct PnetCDF to call MPI collective APIs to
    read and write the file header. The default is "false", meaning the file
    header is only read/written by rank 0, using MPI independent read and write
    APIs. When set to "true", collective APIs will be used with all other
    processes making zero-size requests. See
    [PR #104](https://github.com/Parallel-NetCDF/PnetCDF/pull/104).

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility program
  + ncvalidator - When the file size is larger and smaller than expected, the
    file may still be a valid netCDF file.
    * When larger, this can happen if opening an existing file that contains no
      variable. Deleting a global attribute reduces the file header size. The
      file is still a valid netCDF file.
      [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99) detects
      this mismatch and truncates the file size to the header size.
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
  + File stripe size is no longer used to set the size of file header extent.
    Automatically set the size of file header extent using the file striping
    size obtained from MPI-IO hint 'striping_unit' may yield an unexpectedly
    large file header extent and cause movement of data sections if new
    metadata is added when the program re-enter the define mode.
    See [PR #124](https://github.com/Parallel-NetCDF/PnetCDF/pull/124) and
    [PR #125](https://github.com/Parallel-NetCDF/PnetCDF/pull/125).
  + Use unsigned int to do byte swap.
    See [PR #113](https://github.com/Parallel-NetCDF/PnetCDF/pull/113).
  + Silence Intel icc compilation warnings: when CFLAGS contains
    "-Wimplicit-const-int-float-conversion" and "-Wstringop-overread".
    See [PR #110](https://github.com/Parallel-NetCDF/PnetCDF/pull/110).
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
  + Fix Fortran APIs internal call to nfmpi_inq_buffer_size.
    See [PR #111](https://github.com/Parallel-NetCDF/PnetCDF/pull/111).
  + Fix ncmpi_inq_num_rec_vars and ncmpi_inq_num_fix_vars when opening an
    existing file. See
    [PR #103](https://github.com/Parallel-NetCDF/PnetCDF/pull/103).
  + ncmpidiff -  when checking the dimensions defined in the second files
    whether they are defined in the first file. See 88cd9c1.
  + Use Fortran subroutine Get_Environment_Variable instead of getenv.
    See commit a0b8aca, b796759.

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/tst_redefine.c - test multiple entries of ncmpi__enddef
    [PR #126](https://github.com/Parallel-NetCDF/PnetCDF/pull/126).
  + test/testcases/tst_symlink.c - test NC_CLOBBER on a symbolic link.
  + test/testcases/tst_del_attr.c - test delete attributes. See
    [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99).
  + test/testcases/test_get_varn.c - test get_varn API. See
    [PR #90](https://github.com/Parallel-NetCDF/PnetCDF/pull/90).
  + test/testcases/flexible_var.c - test flexible var API
  + test/testcases/flexible_api.f - test flexible API when bufcount == -1
  + test/testcases/scalar.c - add tests for scalar variables using nonblocking
    APIs. See commit 07ff7b1
  + test/nonblocking/test_bputf.f90, test/nonblocking/test_bputf77.f -
    add tests of APIs inq_buffer_usage and inq_buffer_size. See commit 94ce438

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
  + Hints nc_header_align_size, nc_var_align_size, and nc_record_align_size are
    to align the file header extent, starting file offsets of data section, and
    record variable section, respectively. Note that nc_var_align_size and
    nc_record_align_siz are not to align individual variables. They are
    equivalent to arguments v_align and r_align, respevtively in API
    ncmpi__enddef().
  + Using NC_CLOBBER in ncmpi_create() can be expensive if the file already
    exists. If the existing file is a regular file, then PnetCDF will delete it
    with a call to unlink() first and re-created it later. Calling unlink() on
    may be expensive for some parallel file systems. If the existing file is a
    symbolic link, then PnetCDF will call truncate() or MPI_File_set_size() to
    truncate the file size to zero. Calling truncate() may also be very
    expensive on some file systems, e.g. Lustre. Sporadically a long time spent
    on unlink() and truncate() was observed on Perlmutter.

