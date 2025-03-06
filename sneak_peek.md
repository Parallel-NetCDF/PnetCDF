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

* Update of PnetCDF hints
  + There are three ways in PnetCDF for user to set hints to adjust file header
    extent size:
    1. through setting hints `nc_header_align_size`, `nc_var_align_size`, and
       `nc_record_align_size` in in the environment variable `PNETCDF_HINTS`.
    2. through arguments `h_minfree`, `v_align`, `v_minfree`, and `r_align` of
       `ncmpi__enddef()`.
    3. through setting hints `nc_header_align_size`, `nc_var_align_size`, and
       `nc_record_align_size` in an MPI info object and passing it to calls
       of `ncmpi_create()` and `ncmpi_open()`.
  + Users may set the same hints through one or more of the above methods.
    When a hint is set more than one time, PnetCDF implements the following
    hint precedence.
    * 1st priority: hints set in the environment variable `PNETCDF_HINTS`, e.g.
      `PNETCDF_HINTS="nc_var_align_size=1024"`. Making this the first priority
      is because it allows the same application executable without source code
      modification to run using different settings through a run-time
      environment varaible.
    * 2nd priority: hints used in the arguments of `ncmpi__enddef()`, e.g.
      `ncmpi__enddef(..., v_align=1024,...)`. With the same reason as described
      above, application source codes making a call to APIs using specific
      arguments should have a priority lower than the ones set in a run-time
      environment variable, so that users can try different hints using the
      same executable.
    * 3rd priority: hints set in the MPI info object passed into calls of
      `ncmpi_create()` and `ncmpi_open()`, e.g.
      `MPI_Info_set("nc_var_align_size", "1024");`. This is because MPI info
      can only be passed at the file creation or open time whose hints are
      generally applied to the whole operation to the fille until it is closed.
      However, `ncmpi__endde()` can be called multiple times between a file's
      open and close, as users may want to use a settings unique for individual
      call to `ncmpi__endde()`. Thus `ncmpi__endde()` should have a higher
      priority than MPI info.
  + PnetCDF I/O hint `nc_header_align_size` is essentially the same as hint
    `nc_var_align_size`, but its name is closer to the hint's intent, i.e.
    adjusting the header to reserve space for its growth in the future when new
    data objects are added. Please note when both hints are set, only hint
    `nc_var_align_size` will take effect in PnetCDF and `nc_header_align_size`
    ignored.
  + When there is no fix-sized variable (i.e. non-record variable) defined,
    argument `v_minfree` passsed to `ncmpi__enddef()` is ignored. In this
    case, users should set `h_minfree`, if an extra space is desired.
  + When there is no fix-sized variables defined and none of hints
    `nc_header_align_size`, `nc_var_align_size`, or argument `v_align` is set,
    `nc_record_align_size` or `r_align`, if set, will be used to align the
    header extent.
  + For the above update to the hints used for file layout alignment, see
    PnetCDF See [PR #173](https://github.com/Parallel-NetCDF/PnetCDF/pull/173).

* Clarifications
  + Hint precedence change. PnetCDF version 1.14.0 and prior use a hint
    precedence of `PNETCDF_HINTS` > `MPI info` > `ncmpi__enddef()`. Starting
    from 1.14.1, the precedence has been changed to `PNETCDF_HINTS` >
    `ncmpi__enddef()` > `MPI info`.

