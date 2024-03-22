------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Flexible APIs now can be used as high-level APIs, when argument `bufcount`
    is `NC_COUNT_IGNORE` and `buftype` is an MPI predefined data type. See
    [PR #82](https://github.com/Parallel-NetCDF/PnetCDF/pull/82). Below is an
    example of writing from a memory buffer of type float.
    ```
    ncmpi_put_vara_all(ncid, varid, start, count, buf, NC_COUNT_IGNORE, MPI_FLOAT);
    ```
    is equivalent to
    ```
    ncmpi_put_vara_float_all(ncid, varid, start, count, buf);
    ```
  + PnetCDF now allows a single read/write request from a process of size
    larger than 2 GiB. Such large requests are now passed down to the MP-IO
    library. This change is because MPI 4.0 introduces the large count feature,
    including MPI_Count data type, MPI_XXX_c and MPI_XXX_x APIs that use 8-byte
    integer type to enable large MPI operations. As some MPI libraries today
    have implemented this feature, PnetCDF can now take advantage of it to
    support large single requests. Because of this change, configure option
    `--enable-large-single-req` is thus deprecated. See
    See [PR #131](https://github.com/Parallel-NetCDF/PnetCDF/pull/131)

* New optimization
  + none

* New Limitations
  + Hint `nc_header_read_chunk_size` is limited to `NC_MAX_INT`. PnetCDF reads
    file header in chunks. This hint customizes the chunk size.

* Configure options
  + `--enable-large-single-req` has been removed, as PnetCDF now allows a
    single reqd/write request of size larger than 2 GiB.
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
  + Upgrade `config.guess` and `config.sub` to 2024-01-01.
    See [PR #116](https://github.com/Parallel-NetCDF/PnetCDF/pull/116)
  + FLASH-IO benchmark - add compile flag "-fallow-argument-mismatch" for GNU
    Fortran 10 and later.
    See [PR #114](https://github.com/Parallel-NetCDF/PnetCDF/pull/114)
  + Handle the case when MPICC environment variable is not set and `--with-mpi`
    is not used.  See commit 6142135.
  + Upgrade Autotools version requirement to autoconf 2.71, automake 1.16.5, and
    libtool 2.4.6. (Note this change affects PnetCDF developers only.)
    See [PR #95](https://github.com/Parallel-NetCDF/PnetCDF/pull/95)
    Thanks to Blaise Bourdin for pointing out in
    [Issue #94](https://github.com/Parallel-NetCDF/PnetCDF/issues/94)
    that configure failed when using Intel OneAPI 2022.2.0 compilers. The fix
    is to use autoconf 2.70 and newer.

* New constants
  + `NC_COUNT_IGNORE` - This is used in flexible APIs. When argument `bufcount`
    is `NC_COUNT_IGNORE`, `buftype` must be a predefine MPI datatype and the
    APIs operate as the high-level APIs. Fortran equivalents are
    `NF_COUNT_IGNORE` and `NF90_COUNT_IGNORE`.

* New APIs
  + none

* API syntax changes
  + none

* API semantics updates
  + File open flag `NC_SHARE` is now deprecated. It is still defined, but takes
    no effect.
  + `NC_SHARE` alone is not sufficient to provide data consistency for accessing
    a shared file in parallel and thus is now deprecated.  Because PnetCDF
    follows the MPI file consistency, which only addresses the case when all
    file accesses are relative to a specific file handle created from a
    collective open, `NC_SHARE` becomes invalid. See doc/README.consistency.md
    for more information.

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + none

* New PnetCDF hint
  + `nc_hash_size_dim`:   Set hash table size for dimension names. Default: 256
  + `nc_hash_size_var`:   Set hash table size for variable names. Default: 256
  + `nc_hash_size_gattr`: Set hash table size for global attribute names.
    Default: 64
  + `nc_hash_size_vattr`: Set hash table size for variable attribute names.
    Default: 8
  + The above 4 new hints allow users to set different hash table sizes for
    different objects.  For instance, when the number of variables to be
    defined is large and the number of attributes per variable is small,
    increasing `nc_hash_size_var` can speed up the definition time, and
    reducing `nc_hash_size_vattr` can reduce the memory footprint. See
    [PR #132](https://github.com/Parallel-NetCDF/PnetCDF/pull/132).

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility program
  + `ncvalidator` - When the file size is larger and smaller than expected, the
    file may still be a valid netCDF file.
    * When larger, this can happen if opening an existing file that contains no
      variable. Deleting a global attribute reduces the file header size. The
      file is still a valid netCDF file.
      [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99) detects
      this mismatch and truncates the file size to the header size.
    * When smaller, this can happen if the last variable is partially written.
      The expected file size is calculated based on the full sizes of all
      variables. The file is still a valid netCDF file.
    * In the former case, PR #99 changes `ncvalidator` to report a warning,
      rather than an error.
    * Add printing of the dimension size of a variable when its size is larger
      than the limitation allowed by the file format. See commit 5584d44.
  + Add file src/utils/README.md which gives short descriptions of the utility
    programs and collapsible bullets to display their manual pages.

* Other updates:
  + When file header extent size grows, use 64 MiB per process as the move unit
    size. See [PR #137](https://github.com/Parallel-NetCDF/PnetCDF/pull/137)
  + Since version 1.1.0, PnetCDF has been using file striping size, if
    obtainable from hint `striping_unit` set by users or MPI-IO underneath, to
    align the starting file offset of data section. This offset is also
    referred to as the file header extent, which can be larger than the header
    size to allow header to grow if new data objects are added. Starting from
    this release of 1.13.0, file stripe size is no longer used for this
    purpose. This is because automatically setting file header extent using the
    file striping size may yield an unexpectedly large file header extent and
    cause movement of data sections if new metadata is added.
    See [PR #124](https://github.com/Parallel-NetCDF/PnetCDF/pull/124) and
    [PR #125](https://github.com/Parallel-NetCDF/PnetCDF/pull/125).
  + Use unsigned int to do byte swap.
    See [PR #113](https://github.com/Parallel-NetCDF/PnetCDF/pull/113).
  + Silence Intel icc compilation warnings: when CFLAGS contains
    "-Wimplicit-const-int-float-conversion" and "-Wstringop-overread".
    See [PR #110](https://github.com/Parallel-NetCDF/PnetCDF/pull/110).
  + In all previous PnetCDF's implementations, file header is always
    written/read by rank 0 using MPI independent APIs. This can nullify ROMIO
    hint `romio_no_indep_rw` if set by the user. To warrant no independent
    read/write, PnetCDF now first checks hint `romio_no_indep_rw` and if set to
    `true`, then all file header I/Os are done using MPI collective I/O calls,
    where only rank 0 makes non-zero length requests while all others zero
    length (in order to participate the collective calls). See
    [PR #104](https://github.com/Parallel-NetCDF/PnetCDF/pull/104) and
    [PR #138](https://github.com/Parallel-NetCDF/PnetCDF/pull/138).
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
  + Fix hint values that are actually used. See commit 41e8ef8.
  + Fix residual values of `v_align` and `r_align` when re-entering the define
    mode multiple times.
    See [PR #126](https://github.com/Parallel-NetCDF/PnetCDF/pull/126).
  + Fix Fortran APIs internal call to `nfmpi_inq_buffer_size`.
    See [PR #111](https://github.com/Parallel-NetCDF/PnetCDF/pull/111).
  + Fix `ncmpi_inq_num_rec_vars()` and `ncmpi_inq_num_fix_vars()` when opening
    an existing file. See
    [PR #103](https://github.com/Parallel-NetCDF/PnetCDF/pull/103).
  + `ncmpidiff` -  when checking the dimensions defined in the second files
    whether they are defined in the first file. See 88cd9c1.
  + Use Fortran subroutine `Get_Environment_Variable` instead of `getenv`.
    See commit a0b8aca, b796759.

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/largefile/tst_hash_large_ndims.c - test hashing performance when
    the number of dimensions is large.
  + test/largefile/tst_hash_large_nvars.c - test hashing performance when
    the number of variables is large.
  + test/largefile/tst_hash_large_ngattr.c - test hashing performance when
    the number of global attributes is large.
  + test/largefile/large_header.c - test file header size larger than 2 GiB.
  + test/largefile/large_reqs.c - test a single read/write request of size
    larger than 2 GiB.
  + test/testcases/tst_redefine.c - test multiple entries of `ncmpi__enddef`
    [PR #126](https://github.com/Parallel-NetCDF/PnetCDF/pull/126).
  + test/testcases/tst_symlink.c - test `NC_CLOBBER` on a symbolic link.
  + test/testcases/tst_del_attr.c - test delete attributes. See
    [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99).
  + test/testcases/test_get_varn.c - test `get_varn` API. See
    [PR #90](https://github.com/Parallel-NetCDF/PnetCDF/pull/90).
  + test/testcases/flexible_var.c - test flexible var API
  + test/testcases/flexible_api.f - test flexible API when `bufcount == -1`
  + test/testcases/scalar.c - add tests for scalar variables using nonblocking
    APIs. See commit 07ff7b1
  + test/nonblocking/test_bputf.f90, test/nonblocking/test_bputf77.f -
    add tests of APIs `inq_buffer_usage` and `inq_buffer_size`.
    See commit 94ce438

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
  + Hints `nc_header_align_size` and `nc_record_align_size` are to align the
    file header extent and starting file offset of the record variable section,
    respectively. Hint `nc_record_align_size` is not to align the offsets of
    individual record variables.
  + Prior to version 1.13.0, hint `nc_var_align_size` is used to align the
    starting file offsets of individual fixed-size variables. Beginning from
    version 1.13.0, this hint is used to align only the starting file offset of
    the entire data section, not individual variables. The reason of such
    change is because PnetCDF nonblocking APIs can aggregate multiple I/O
    requests and MPI-IO today has been optimized to align I/O to file striping
    boundaries, which makes aligning starting offsets of individual variables
    less effective and may create empty space in the file.
  + Using `NC_CLOBBER` in `ncmpi_create()` can be expensive if the file already
    exists. If the existing file is a regular file, then PnetCDF will delete it
    with a call to `unlink()` first and re-created it later. Calling `unlink()` on
    may be expensive for some parallel file systems. If the existing file is a
    symbolic link, then PnetCDF will call `truncate()` or `MPI_File_set_size()` to
    truncate the file size to zero. Calling `truncate()` may also be very
    expensive on some file systems, e.g. Lustre. Sporadically a long time spent
    on `unlink()` and `truncate()` was observed on Perlmutter.

