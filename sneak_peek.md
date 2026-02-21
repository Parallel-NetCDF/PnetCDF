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

* Configure updates:
  + none

* New constants
  + none

* New APIs
  + none

* API syntax changes
  + none

* API semantics updates
  + API ncmpi_inq_header_size() now can be called in the define mode. This API
    returns the file header size with metadata defined by the time of the call.
    This information can be helpful to pick proper values for arguments
    h_minfree, v_align, v_minfree, r_align when calling ncmpi__enddef() to
    allocate a sufficiently large free space for file header extent and
    variable data sections to grow without moving data already stored in the
    file, i.e. when adding new variables, dimensions, or attributes.
    See [PR #201](https://github.com/Parallel-NetCDF/PnetCDF/pull/201).

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + none

* New PnetCDF hints
  + 'nc_data_move_chunk_size' -- When adding new data objects into an existing
    file, the data section may need to be moved to a higher file offset. The
    movement is performed in chunks. This hint allows users to customized the
    chunk size. The default is 1048576 bytes, i.e. 1 MiB.
    See [PR #203](https://github.com/Parallel-NetCDF/PnetCDF/pull/203).
  + 'nc_striping' -- When creating a new file on the Lustre file system, this
    hint advises PnetCDF to set the file's striping configuration. The hint
    value is either "auto" or "inherit". The former sets the new file's
    striping unit to 1MiB and striping count to the number of compute nodes
    found in the MPI communicator passed to `ncmpi_create()`. The latter sets
    the stripings of new file to inherit the folder's striping settings, if the
    folder's striping is set. However, if users also set hint `striping_factor`
    or `striping_unit`, then their values will be used instead of parent
    folder's. This hint's default value is "auto".
    See [PR #222](https://github.com/Parallel-NetCDF/PnetCDF/pull/222).

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + none

* Other updates:
  + none

* Bug fixes
  + Fix data movement when new record variables are added to an existing file
    that does not change the starting offset of record variable section.
    See [PR #199](https://github.com/Parallel-NetCDF/PnetCDF/pull/199).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/tst_grow_data.c -- adding new variables by re-entering the
    define mode multiple time, but does not cause file header extent to grow.
    It also tests a case when adding a new record variable that does not change
    the starting offset of the record variable section in the file.

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
  + none

