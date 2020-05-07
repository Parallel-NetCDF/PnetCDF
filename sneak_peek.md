------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

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
  + ncvalidator now reports the name of variable that violates the NetCDF
    limitation on large variables for CDF-2 files
  + add corrupted file bad_large_fixed_var.nc2 that contains one large
    fixed-size variables that is not defined last
  + add corrupted file bad_large_rec_2_vars.nc2 that contains 2 large record
    variables
  + add corrupted file bad_large_rec_var.nc2 that contains 1 large record
    variable that is not defined last
  + add URLs of NetCDF limitations on large variables for CDF-1 and 2 file
    formats

* Other updates:
  + When using NC_CLOBBER in the call of ncmpi_create(), use access() to check
    whether file exists. If not, successive calls to truncate() or unlink() can
    be skip.
  + Improve detection of HDF5 signature. The HDF5 signature is located at the
    beginning of the HDF5 superblock, but the location of HDF5 superblock may
    not be at the beginning of the file. It is located at byte offset 0, byte
    offset 512, and at successive locations in the file, each a multiple of two
    of the previous location; in other words, at these byte offsets: 0, 512,
    1024, 2048, and so on.

* Bug fixes
  + Fix NC_CLOBBER mode for ncmpi_create() when called on existing symbolically
    linked files. In all previous PnetCDF implementations, symbolic links, like
    other regular files, was first deleted by unlink() and then created. This
    can result in an unexpected outcome, i.e. the symbolic link being deleted.
    NetCDF library implements this differently, by adding O_TRUNC flag when
    calling open() to truncate the file to length 0. Historically, PnetCDF did
    not adopt the same approach because MPI does not define a similar flag to
    O_TRUNC and the only way to achieve the file clobber effect is to through
    MPI_File_set_size(), which can be expensive as the function takes an MPI
    file handler argument, which requires to open the file first with a call to
    MPI_File_open().
  + Fix various compile and link bugs when NAG Fortran is used. Bugs include
    flag needed to verbose linking output, unrecognized link option -pthread,
    unmatched C compiler underneath. Thanks Sergey Kosukhin for providing the
    fix in [PR #59](https://github.com/Parallel-NetCDF/PnetCDF/pull/59)
    and [PR #60](https://github.com/Parallel-NetCDF/PnetCDF/pull/60)
  + Fix a bug of calling Fortran getarg() with the first argument k with a
    value > 0 when there is no command-line argument is used. NAG Fortran may
    crash the program. See
    [f16bd3c](https://github.com/Parallel-NetCDF/PnetCDF/commit/f16bd3c1ba1b08eade2384f094c519f3f2dc114e)
  + Fix a bug that limits FLASH-IO to run on less than 16K MPI processes. See
    [1d84fa5](https://github.com/Parallel-NetCDF/PnetCDF/commit/1d84fa5d54ca9179da4a5b1a4ee3b92cc92287ed)

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
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

