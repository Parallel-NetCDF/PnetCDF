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
  + none

* Bug fixes
  + Fix various compile and link bugs when NAG Fortran is used. Bugs include
    flag needed to verbose linking output, unrecognized link option -pthread,
    unmatched C compiler underneath. Thanks Sergey Kosukhin for providing the
    fix in [PR #59](https://github.com/Parallel-NetCDF/PnetCDF/pull/59)
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

