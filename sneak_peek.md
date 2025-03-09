------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + none

* New optimization
  + When file header extent size grows, moving the data section to a higher
    file offset has changed to be done in chunks of 16 MB per process.
    See [PR #174](https://github.com/Parallel-NetCDF/PnetCDF/pull/174),

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

* Clarifications about of PnetCDF hints
  + There are three ways in PnetCDF for user to set hints to align the starting
    file offset for the data section (header extent) and record variable
    section.
    1. through a call to API `nc_header_align_size` by setting arguments of
       `h_minfree`, `v_align`, `v_minfree`, and `r_align`.
    2. through an MPI info object passed to calls of `ncmpi_create()` and
       `ncmpi_open()`. Hints are `nc_header_align_size`, `nc_var_align_size`,
       and `nc_record_align_size`.
    3. through a run-time environment variable `PNETCDF_HINTS`. Hints are
       `nc_header_align_size`, `nc_var_align_size`, and `nc_record_align_size`.
  + As the same hints may be set by one or more of the above methods, PnetCDF
    implements the following hint precedence.
    * `PNETCDF_HINTS` > `ncmpi__enddef()` > `MPI info`.
    * 1st priority: hints set in the environment variable `PNETCDF_HINTS`, e.g.
      `PNETCDF_HINTS="nc_var_align_size=1048576"`. Making this the first
      priority is because it allows to run the same application executable
      without source code modification using different alignment settings
      through a run-time environment variable.
    * 2nd priority: hints set in the MPI info object passed to calls of
      `ncmpi_create()` and `ncmpi_open()`, e.g.
      `MPI_Info_set("nc_var_align_size", "1048576");`. The reasoning is when a
      3rd-party library built on top of PnetCDF implements its codes using
      'ncmpi__enddef'. An application that uses such 3rd-party library can pass
      an MPI info object to it, which further passes the info to PnetCDF. This
      precedence allows that application to exercise different hints without
      changing the 3rd-party library's source codes.
    * 3rd priority: hints used in the arguments of `ncmpi__enddef()`, e.g.
      `ncmpi__enddef(..., v_align=1048576,...)`.
  + PnetCDF I/O hint `nc_header_align_size` is essentially the same as hint
    `nc_var_align_size`, but its name appears to be closer to the hint's
    intent, i.e. to reserve some space for the header growth in the future when
    new data objects are added. Please note when both hints are set, only hint
    `nc_var_align_size` will take effect and `nc_header_align_size` ignored.
  + When there is no fix-sized variable (i.e. non-record variable) defined,
    argument `v_minfree` passed to `ncmpi__enddef()` is ignored. In this
    case, users should set `h_minfree`, if an extra header space is desired.
  + When there is no fix-sized variables defined and none of hints
    `nc_header_align_size`, `nc_var_align_size`, or argument `v_align` is set,
    `nc_record_align_size` or `r_align`, if set, will be used to align the
    header extent.
  + For the above update to the hint precedence, see
    PnetCDF See [PR #173](https://github.com/Parallel-NetCDF/PnetCDF/pull/173).

