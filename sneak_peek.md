------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New feature
  + Intra-node aggregation for read requests is added. This is the complement
    to the write counterpart first added to version 1.14.0. Now the intra-node
    aggregation feature supports both write and read operations. This feature
    can be enabled by setting hint `nc_num_aggrs_per_node` to the desired
    number of aggregators per compute node.

* New optimization
  + When creating a new file on the Lustre file system, PnetCDF will try to set
    the Lustre file striping count to the number of compute nodes (NUMA nodes)
    allocated to the MPI communicator passed to the call of `ncmpi_create()`.
    See the section of new PnetCDF hints below for other user options to set
    the new file's striping configuration.

* New I/O driver
  + A submodule named "GIO" has been added as the PnetCDF's default I/O driver.
    GIO is a parallel I/O library, sharing the design concept of the MPI-IO.
    GIO implements several performance improvement, such as:
    * automatically sets `cb_nodes` based on the number of compute nodes and
      number of MPI processes running per node.
    * supports hint `cb_buffer_size` to be a multiple of file striping unit
      size for Lustre.
    * supports Lustre overstriping through hint `overstriping_ratio`. See the
      "New PnetCDF hints" section below.

* New Limitations
  + none

* Configure options

* Configure updates:
  + none

* New constants
  + none

* New APIs
  + none

* APIs deprecated
  + The "vard" family APIs introduced in version 1.6.0 are now deprecated.
    These are the API family that take an argument of MPI derived data type
    describing the file access layout, which is used as the fileview by the
    underlying MPI-IO library. The deprecation is because direct file access
    can be a security risk and error prone.
    * ncmpi_get_vard()
    * ncmpi_get_vard_all()
    * ncmpi_put_vard()
    * ncmpi_put_vard_all()

* API syntax changes
  + none

* API semantics updates
  + API `ncmpi_inq_header_size()` now can be called in the define mode. This
    API returns the file header size with metadata defined by the time of the
    call. The inquired file header size can be helpful to pick proper values
    for arguments `h_minfree`, `v_align`, `v_minfree`, `r_align` when calling
    API `ncmpi__enddef()` which can allow to preserve a sufficiently large free
    space for file header extent and variable data sections to grow later
    without moving data already stored in the file, i.e. when adding new
    variables, dimensions, or attributes.
    See [PR #201](https://github.com/Parallel-NetCDF/PnetCDF/pull/201).

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + `NC_EFSTYPE` indicates an error when an invalid file system type is
    detected by the GIO library.
  + `NC_EDRIVER` indicates an invalid PnetCDF I/O driver is set in I/O hint.
    The current valid drivers are "gio" and "mpiio".
  + `NC_EFILEVIEW` indicates a PnetCDF internal error when creating an MPI
    fileview whose offsets violate the MPI standard requirement of being in a
    monotonically non-decreasing order.

* New PnetCDF hints
  + `nc_data_move_chunk_size` -- When adding new data objects into an existing
    file, the data sections may need to be moved to a higher file offset. The
    movement is performed in chunks. This hint allows users to customized the
    chunk size. The default is 1048576 bytes, i.e. 1 MiB.
    See [PR #203](https://github.com/Parallel-NetCDF/PnetCDF/pull/203).
  + `file_striping` -- When creating a new file on the Lustre file system, this
    hint advises PnetCDF to set the new file's striping configuration. The hint
    value is either "auto" or "inherit". The default is "auto". The former sets
    the new file's striping unit size to 1 MiB and striping count to the number
    of compute nodes found in the MPI communicator passed to the call of
    `ncmpi_create()`.  The latter sets the striping to inherit the new file's
    parent folder's settings, if the folder's striping has been set. Note that
    if users also set the MPI-IO hint `striping_factor` or `striping_unit`,
    then these MPI-IO hints will take a higher precedence.
    See [PR #222](https://github.com/Parallel-NetCDF/PnetCDF/pull/222).
  + `nc_driver` -- To select an I/O driver. A string value of "gio" is to use
    the GIO driver, an external library built as a sub-module of PnetCDF, and
    "mpiio" the MPI-IO driver. The default is "gio".
  + `overstriping_ratio` -- This hint instructs GIO driver to make use of
    Lustre's overstriping feature which allows to create a new file to have
    more than one stripe per OST. For instance, when hints `striping_factor`
    and `overstriping_ratio` are set to 16 and 2, respectively, the newly
    created file will be striped across 16/2=8 OSTs.

* New run-time environment variables
  + `GIO_ONLY` - when testing PnetCDF by running command `make check` and `make
    ptest`, setting this environment variable will run all the tests using the
    default GIO driver and skip MPI-IO driver. Without setting this environment
    variable, the default tests both I/O drivers.

* Build recipes
  + `doc/NERSC.md` - note for building PnetCDF on machines at NERSC

* Updated utility programs
  + `pnetcdf_version` - add information about MPI compiler vendors, supported
    MPI standard version, names of base compilers, and LMOD PrgEnv module
    loaded if available.
  + `ncmpidiff` now will error out when the two input file names are identical.

* Other updates:
  + none

* Bug fixes
  + Fix problem of building PnetCDF shared library when using NVIDIA HPC
    compilers, nvc, nvc++, and nvfortran.
  + Fix data movement when new record variables are added to an existing file
    that does not change the starting offset of record variable section.
    See [PR #199](https://github.com/Parallel-NetCDF/PnetCDF/pull/199).

* New example programs
  + none

* New I/O benchmarks
  + benchmarks/C/indep_data_obj_create.c - It evaluates the performance of
    creating a large number of data objects (dimensions and variables) in an
    output file. This benchmark is used in a paper published in the PDSW
    Workshop of SC 2025. DOI: 10.1145/3731599.3767512

* New test programs
  + test/testcases/tst_grow_data.c -- It tests the case when adding new
    variables by re-entering the define mode multiple time without causing the
    file header extent to grow. It also tests a case when adding a new record
    variable that does not change the starting offset of the record variable
    section in the file.

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

