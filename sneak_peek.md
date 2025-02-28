------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + none

* New Limitations
  + none

* Configure options
  + For PnetCDF developers, the requirement for libtool version has been
    changed to 2.5.4, due to an issue on Mac OS when using OpenMPI. See
    [Issue #155](https://github.com/Parallel-NetCDF/PnetCDF/issues/155),
    [Issue #163](https://github.com/Parallel-NetCDF/PnetCDF/issues/163),
    and [PR #164](https://github.com/Parallel-NetCDF/PnetCDF/pull/164).

* Configure updates:
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

* New PnetCDF hints
  + none

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + none

* Other updates:
  + none

* Bug fixes
  + Fix setting of user hint nc_ibuf_size.
    See [PR #161](https://github.com/Parallel-NetCDF/PnetCDF/pull/161).

* New example programs
  + none

* New programs for I/O benchmarks
  + WRF-IO contains an extraction of the I/O kernel of WRF (Wether Research
    and Forecast Model, a weather prediction computer simulation program
    developed at NCAR) that can be used to evaluate the file write performance
    of WRF. It's data partitioning pattern is a 2D block-block checkerboard
    pattern, along the longitude and latitude.
    See [PR #165](https://github.com/Parallel-NetCDF/PnetCDF/pull/165).

* New test program
  + none

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
  + There are three ways in PnetCDF for user to set hints to adjust file header
    extent size, for example. Users may set the same hints multiple times
    during the run of their applications. When this happens, PnetCDF implements
    the following hint precedence.
    * 1st priority: hints set in the environment variable `PNETCDF_HINTS`, e.g.
      `PNETCDF_HINTS="nc_var_align_size=1024"`
    * 2nd priority: hints passed from arguments of `ncmpi__enddef()`, e.g.
      `ncmpi__enddef(..., v_align=1024,...)`
    * 3rd priority: hints set in the MPI info objects passed into calls to
      `ncmpi_create()` and `ncmpi_open()`, e.g.
      `MPI_Info_set("nc_var_align_size", "1024");`
    See [PR #173](https://github.com/Parallel-NetCDF/PnetCDF/pull/173).
  + PnetCDF I/O hint `nc_header_align_size` is essentially the same as hint
    `nc_var_align_size`, but name of the former is closer to the intent, i.e.
    to adjust the header space to reserve space for possible expansion in the
    future when new data objects are added. However, when both hints were set
    by the users, only hint `nc_var_align_size` will take effect in PnetCDF and
    `nc_header_align_size` will be ignored.

