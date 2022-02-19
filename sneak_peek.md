------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + Improve the performance when nonblocking requests contain more than one
    record.
    See [a00f5d2](https://github.com/Parallel-NetCDF/PnetCDF/commit/a00f5d2).

* New Limitations
  + none

* Update configure options
  + Retire configure options `--enable-netcdf4` and `--enable-adios`, as there
    are already options `--with-netcdf4` and `--with-adios`. According to the
    autoconf manual, `--enable-feature` is for internal packages and
    `--with-feature` is for external.  Adding `--with-feature` is equivalent to
    enabling the feature.

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
  + Replace the definition of "difference ratio" used in utility programs
    `cdfdiff` and `ncmpidiff` with formula
    ```
    |x - y| / max(|x|, |y|)
    ```
    where |x| means the absolute value of x.
    See [issue #78](https://github.com/Parallel-NetCDF/PnetCDF/issues/78) and
    [9b165ceb](https://github.com/Parallel-NetCDF/PnetCDF/commit/9b165ceb).
  + Utility programs `cdfdiff` and `ncmpidiff` now report only the first
    variable element that are different or fail to meet the tolerances. Their
    values and differences are now printed on the standard output.
    See [75d20e6](https://github.com/Parallel-NetCDF/PnetCDF/commit/75d20e6)
    and [58d2d17](https://github.com/Parallel-NetCDF/PnetCDF/commit/58d2d17).
  + Use the same user-provided tolerant difference and tolerant difference
    ratio (through command-line option '-t') to check all variables.
    See [issue #78](https://github.com/Parallel-NetCDF/PnetCDF/issues/78) and
    [75d20e6](https://github.com/Parallel-NetCDF/PnetCDF/commit/75d20e6).

* Other updates:
  + IBM XLF compiler on Summit at OLCF requires compile flag "-qfixed" when
    compiling Fortran programs written in fixed form. Thus, the fixed form flag
    detected at configure time has been added to `AM_FFLAGS` when compiling
    fixed form programs. Similarly, the free form flag has been added to
    `AM_FCFLAGS` when compiling free form programs. This issue affects only the
    test and example Fortran programs. The PnetCDF library is intact.
    See [PR #73](https://github.com/Parallel-NetCDF/PnetCDF/pull/73)
  + Add all PnetCDF I/O hints to the inquired MPI info object returned by API
    `ncmpi_inq_file_info()`.
    See [f0e65cf](https://github.com/Parallel-NetCDF/PnetCDF/commit/f0e65cf).

* Bug fixes
  + Fix a bug in utility program `ncvalidator` for the case when there is no
    record written, but both fixed-size and record variables are defined.
    See [b34bfcd](https://github.com/Parallel-NetCDF/PnetCDF/commit/b34bfcd).
  + Fix configure bug of setting environment variable `SEQ_CC` to the
    sequential `CC` extracted from `MPICC`. Also add configure help message for
    environment variable `SEQ_CC`. Thanks to Carl Ponder for reporting.
    See [4978f6d](https://github.com/Parallel-NetCDF/PnetCDF/commit/4978f6d).
  + Fix a bug when C function `truncate` is not available. Argument `fh` of
    `MPI_File_close` is a pointer.
    See [3e331a6](https://github.com/Parallel-NetCDF/PnetCDF/commit/3e331a6).
  + When using an MPI compiler whose Fortran feature was disabled, the MPI
    Fortran constants and datatypes may not be defined in the header file
    `mpi.h`. This is the case for Open MPI (tested with 4.0.2). PnetCDF used
    some Fortran datatypes without checking whether they are defined, which can
    fail when running 'make'. A fix has been added to check whether the Fortran
    feature is disabled and wraps around the Fortran datatypes with C directive
    'ifdef ENABLE_FORTRAN'. Thanks Bert Wesarg for reporting this bug.
    See [issue #68](https://github.com/Parallel-NetCDF/PnetCDF/issues/68) and
    [PR #69](https://github.com/Parallel-NetCDF/PnetCDF/pull/69).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + none

* Issues with NetCDF library
  + Test program [test/nc4/tst_rec_vars.c](test/nc4/tst_rec_vars.c) fails to
    run when using NetCDF 4.8.0 and 4.8.1. Thanks Bruno Pagani for reporting.
    See [issue #72](https://github.com/Parallel-NetCDF/PnetCDF/issues/72).

* Conformity with NetCDF library
  + none

* Discrepancy from NetCDF library
  + none

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications
  + Nonblocking APIs have yet been updated to support subfiling feature.

