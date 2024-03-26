------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + A single read/write request made by an MPI process is now allowed to be of
    size larger than 2 GiB. Such large requests will be passed to the MP-IO
    library. This feature makes use of "the large count feature" introduced in
    MPI standard 4.0, which includes `MPI_XXX_c` APIs whose arguments are of
    type `MPI_Count`. `MPI_Count` can be an 8-byte integer type, enabling large
    MPI operations. As some MPI libraries today have begun implementing MPI
    4.0, PnetCDF now can rely on the MPI libraries to support large single
    requests. When the MPI library used to build PnetCDF does not support large
    requests, the MPI errors are returned. Because of this change, the PnetCDF
    configure option `--enable-large-single-req` is thus deprecated.
    See [PR #131](https://github.com/Parallel-NetCDF/PnetCDF/pull/131)
  + Flexible APIs now can operate as high-level APIs, when argument `bufcount`
    is set to `NC_COUNT_IGNORE` and `buftype` is set to an MPI predefined data
    type. See [PR #82](https://github.com/Parallel-NetCDF/PnetCDF/pull/82).
    Below is a write example from a buffer of type float.
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
  + Hint `nc_header_read_chunk_size`, introduced in version 1.4.0, is now
    limited to `NC_MAX_INT`. As PnetCDF reads file header in chunks, this hint
    can be used to customize the chunk size. The default is 256 KB.
    See [4209056](https://github.com/Parallel-NetCDF/PnetCDF/commit/4209056e9a66465421f7ce9f1b44518923638b04)

* Configure options
  + `--enable-large-single-req` is deprecated and removed, as PnetCDF now
    allows a single reqd/write request of size larger than 2 GiB.
  + `--disable-file-sync` is deprecated and removed. This configure option
    alone was not able to provide a sufficient data consistency. Users are
    suggested to call `ncmpi_sync` and `MPI_Barrier` to achieve a desired
    consistency, as suggested by MPI standard.
  + A new option `--enable-install-examples` installs the example programs
    under folder `${prefix}/pnetcdf_examples` along with run script files. An
    example is `${prefix}/pnetcdf_examples/C/run_c_examples.sh`. The default of
    this option is `disabled`.
    See [PR #91](https://github.com/Parallel-NetCDF/PnetCDF/pull/91)
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
    is not used.  See commit
    [6142135](https://github.com/Parallel-NetCDF/PnetCDF/commit/61421356ecd38878a4ef46771ed6520d4257251f)
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
    See [PR #92](https://github.com/Parallel-NetCDF/PnetCDF/pull/92)

* New APIs
  + none

* API syntax changes
  + none

* API semantics updates
  + File open flag `NC_SHARE` is now deprecated. It is still defined, but takes
    no effect.
    See [PR #119](https://github.com/Parallel-NetCDF/PnetCDF/pull/119)
  + `NC_SHARE` alone is not sufficient to provide data consistency for accessing
    a shared file in parallel and thus is now deprecated. PnetCDF follows the
    file consistency defined in MPI standard, which only addresses the case
    when all file accesses are relative to a specific file handle created from
    a collective open, and thus `NC_SHARE` becomes invalid. See
    doc/README.consistency.md for more information.

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + none

* New PnetCDF hints
  + `nc_hash_size_dim` sets the hash table size for dimension names.
    Default: 256
  + `nc_hash_size_var` sets the hash table size for variable names.
    Default: 256
  + `nc_hash_size_gattr` sets the hash table size for global attribute names.
    Default: 64
  + `nc_hash_size_vattr` sets the hash table size for variable attribute names.
    Default: 8
  + The above 4 new hints can be used to set different hash table sizes for
    dimensions, variables, and attributes. Hashing tables are used for quick
    data object name lookup. It can be useful for files containing a large
    number of dimensions, variables, and attributes. For instance, when the
    number of variables to be defined is large and the number of attributes per
    variable is small, increasing the value of `nc_hash_size_var` can speed up
    the variable definition and inquiring time. On the other hand, setting a
    smaller value for hint `nc_hash_size_vattr` can reduce memory footprint.
    See [PR #132](https://github.com/Parallel-NetCDF/PnetCDF/pull/132).

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + `ncvalidator` - When the file size is larger or smaller than what is
    calculated based on the metadata stored in the file header, the file may
    still be a valid netCDF file.
    * The larger-than-expected case can happen if opening an existing file that
      contains no variable. Deleting a global attribute already defined in the
      file will reduce the file header size. In this case, the file is still a
      valid netCDF file. Such mismatch will be detected and the file size
      truncated to the header size.
      See [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99)
    * The smaller-than-expected case can happen if the last variable is
      partially written. The expected file size is calculated based on the full
      sizes of all variables. In this case, the file is still a valid netCDF
      file, and `ncvalidator` will report a warning, rather than an error.
    * Print the dimension size of a variable on stdout when its size is larger
      than the limitation allowed by the file format. See commit
      [5584d44](https://github.com/Parallel-NetCDF/PnetCDF/commit/5584d44a433a68966b0be601e7a73e939c695dbf)
  + Add file src/utils/README.md which gives short descriptions of the utility
    programs and collapsible bullets to display their manual pages.

* Other updates:
  + When file header extent size grows, PnetCDF now uses a movement unit per
    process of size up to 64 MiB.
    See [PR #137](https://github.com/Parallel-NetCDF/PnetCDF/pull/137)
  + Since version 1.1.0, PnetCDF has been using file striping size, if
    obtainable from the MPI-IO hint `striping_unit`, to align the starting file
    offset of the data section. This offset is also referred to as the file
    header extent, which can be larger than the header size to allow header to
    grow when new data objects are added. Starting from this release, file
    stripe size is no longer used for setting the starting offset of the data
    section. This is because automatically setting file header extent using the
    file striping size may grow the file header unexpectedly when adding new
    objects to an existing file.
    See [PR #124](https://github.com/Parallel-NetCDF/PnetCDF/pull/124) and
    [PR #125](https://github.com/Parallel-NetCDF/PnetCDF/pull/125).
  + Use unsigned int to perform byte swap.
    See [PR #113](https://github.com/Parallel-NetCDF/PnetCDF/pull/113).
  + Silence Intel icc compilation warnings: when CFLAGS contains
    "-Wimplicit-const-int-float-conversion" and "-Wstringop-overread".
    See [PR #110](https://github.com/Parallel-NetCDF/PnetCDF/pull/110).
  + In all previous PnetCDF's implementations, file header is always written/
    read by rank 0 using MPI independent APIs. This can nullify ROMIO hint
    `romio_no_indep_rw` if set by the user. To warrant no independent read/
    write user hint, PnetCDF now checks hint `romio_no_indep_rw` and if set to
    `true`, then all file header I/Os are made through MPI collective I/O
    calls, where only rank 0 makes non-zero length requests while all others
    zero length (in order to participate the collective calls). See
    [PR #104](https://github.com/Parallel-NetCDF/PnetCDF/pull/104) and
    [PR #138](https://github.com/Parallel-NetCDF/PnetCDF/pull/138).
  + In all prior versions, the file name was checked whether it contains
    character ':'. The prefix name ending with ':' is considered by ROMIO as
    the file system type name. The prefix name, if found, is then stripped, so
    the file name can be used in the POSIX function calls internally. However,
    the prefix was not checked against the file system type names recognized
    by ROMIO. Starting from this release, the prefix is checked against the
    known file system type names to ROMIO. If the prefix is not one of the
    recognized types, e.g.  "ufs", "nfs", "xfs", "pvfs2", "gpfs", "panfs",
    "lustre", "daos", "testfs", "ime", or "quobyte", then the prefix name is
    not stripped. This change is for the case when the file name contains ':',
    but it is not for specifying the file system type. See
    [PR #79](https://github.com/Parallel-NetCDF/PnetCDF/pull/79) and
    [MPICH PR 5951](https://github.com/pmodels/mpich/pull/5951).

* Bug fixes
  + Return I/O hints that are actually used. See commit
    [0dcf628](https://github.com/Parallel-NetCDF/PnetCDF/commit/0dcf6284106e42faa5e8fb9ab1aa5d52917ff892).
  + Fix residual values of `v_align` and `r_align` when re-entering the define
    mode multiple times.
    See [PR #126](https://github.com/Parallel-NetCDF/PnetCDF/pull/126).
  + Fix Fortran API `nf90mpi_Inq_buffer_size` which should call
    `nfmpi_inq_buffer_size` internally.
    See [PR #111](https://github.com/Parallel-NetCDF/PnetCDF/pull/111).
  + Fix `ncmpi_inq_num_rec_vars()` and `ncmpi_inq_num_fix_vars()` when opening
    an existing file. See
    [PR #103](https://github.com/Parallel-NetCDF/PnetCDF/pull/103).
  + `ncmpidiff` -  when checking the dimensions defined in the second files
    whether are also defined in the first file. See commit
    [88cd9c1](https://github.com/Parallel-NetCDF/PnetCDF/commit/88cd9c187b9b3dc9018b066e905a10d0c74488f8).
  + Use Fortran subroutine `Get_Environment_Variable` instead of `getenv` if it
    is available.  See commits
    [a0b8aca](https://github.com/Parallel-NetCDF/PnetCDF/commit/a0b8acabc3bf7a7a35d878c9db2afb71942bb7a9), and
    [b796759](https://github.com/Parallel-NetCDF/PnetCDF/commit/b796759f5a8e1749e7c168e05fe1665a35e2a2a1).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/largefile/large_attr.c - tests attributes of size > 2 GiB.
  + test/largefile/tst_hash_large_ndims.c - tests hashing performance when
    the number of dimensions is large.
  + test/largefile/tst_hash_large_nvars.c - tests hashing performance when
    the number of variables is large.
  + test/largefile/tst_hash_large_ngattr.c - tests hashing performance when
    the number of global attributes is large.
  + test/largefile/large_header.c - tests file header size larger than 2 GiB.
  + test/largefile/large_reqs.c - tests a single read/write request of size
    larger than 2 GiB.
  + test/testcases/tst_redefine.c - tests multiple entries of `ncmpi__enddef`
    [PR #126](https://github.com/Parallel-NetCDF/PnetCDF/pull/126).
  + test/testcases/tst_symlink.c - tests `NC_CLOBBER` on a symbolic link.
  + test/testcases/tst_del_attr.c - tests delete attributes. See
    [PR #99](https://github.com/Parallel-NetCDF/PnetCDF/pull/99).
  + test/testcases/test_get_varn.c - tests `get_varn` API. See
    [PR #90](https://github.com/Parallel-NetCDF/PnetCDF/pull/90).
  + test/testcases/flexible_var.c - tests flexible var API
  + test/testcases/flexible_api.f - tests flexible API when `bufcount == -1`
  + test/testcases/scalar.c - adds tests for scalar variables using nonblocking
    APIs. See commit
    [07ff7b1](https://github.com/search?q=repo%3AParallel-NetCDF%2FPnetCDF+07ff7b1&type=commits)
  + test/nonblocking/test_bputf.f90, test/nonblocking/test_bputf77.f -
    add tests of APIs `inq_buffer_usage` and `inq_buffer_size`. See commit
    [94ce438](https://github.com/Parallel-NetCDF/PnetCDF/commit/94ce438262fe7fcade031dcae1a677a827549bb3)

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
    respectively.
  + Hint `nc_record_align_size` is not to align the offsets of individual
    record variables.
  + Prior to version 1.13.0, hint `nc_var_align_size` is used to align the
    starting file offsets of individual fixed-size variables. This design was
    to reduce file lock contention in MPI collective I/O operations. Beginning
    from this release (version 1.13.0), this hint is used to align only the
    starting file offset of the entire data section, not individual variables.
    The reason of such change is because PnetCDF nonblocking APIs can aggregate
    multiple I/O requests and MPI-IO today has been optimized to align I/O to
    file striping boundaries, which makes aligning starting offsets of
    individual variables less effective and may create large empty space in the
    file.
  + Using `NC_CLOBBER` in `ncmpi_create()` can be expensive if the file already
    exists. When the existing file is a regular file, PnetCDF will delete it
    with a call to `unlink()` first and re-created it. Calling `unlink()` may
    be expensive for some parallel file systems. When the existing file is a
    symbolic link, PnetCDF will call `truncate()` or `MPI_File_set_size()` to
    truncate the file size to zero. Calling `truncate()` may also be very
    expensive on some file systems, e.g. Lustre. Sporadically a long time spent
    on `unlink()` and `truncate()` was observed on Perlmutter at NERSC.

