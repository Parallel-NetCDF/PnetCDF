------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Intra-node aggregation for write requests -- When the number of MPI
    processes allocated to a compute node is large, this feature selects a
    subset of processes per node to be an I/O aggregator. The non-aggregators
    send their requests to the assigned aggregators, and then the aggregators
    make aggregated requests to the file. This feature can effectively reduce
    communication congestion due to too many pending asynchronous messages
    produced in the collective write inside of MPI-IO. This new feature can be
    enabled by setting the PnetCDF I/O hint 'nc_num_aggrs_per_node' to the
    desired number of aggregators per compute node.
    See [PR #156](https://github.com/Parallel-NetCDF/PnetCDF/pull/156).
  + Support MPI derived data types that are constructed from the large-count
    derived datatype constructors introduced in MPI 4.0.
    See [PR #145](https://github.com/Parallel-NetCDF/PnetCDF/pull/145).

* New optimization
  + When running sequentially (i.e. number of processes is 1), PnetCDF calls
    the MPI independent I/O functions and avoids calls to MPI_Barrier,
    MPI_Bcast, and MPI_Allreduce.
    See [PR #149](https://github.com/Parallel-NetCDF/PnetCDF/pull/149).

* New Limitations
  + none

* Configure options
  + The default has been changed to build both shared and static libraries.
    See [PR #143](https://github.com/Parallel-NetCDF/PnetCDF/pull/143).

* Configure updates:
  + Fix `pnetcdf-config` of reflecting the installation path when installation
    is done by running command `make install DESTDIR=/alternate/directory`
    which prepends '/alternate/directory' before all installation names.
    See [PR #154](https://github.com/Parallel-NetCDF/PnetCDF/pull/154).

* New constants
  + A new C macro `NC_FillValue` replaces `_FillValue` and thus `_FillValue` is
    now deprecated This conforms with NetCDF4's change in its version 4.9.3
    release. See [PR #153](https://github.com/Parallel-NetCDF/PnetCDF/pull/153).

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
  + 'nc_num_aggrs_per_node' -- To enable the intra-node aggregation, this I/O
    hint can set to a positive integral value, which indicates the desired
    number of processes per compute node to be selected as the aggregators.
    Setting it to 0 disables the aggregation, which is also the default mode.
    See [PR #156](https://github.com/Parallel-NetCDF/PnetCDF/pull/156).

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + none

* Other updates:
  + More document for comparing PnetCDF and NetCDF4 has been added to file
    doc/netcdf4_vs_pnetcdf.md
    See [PR #152](https://github.com/Parallel-NetCDF/PnetCDF/pull/152) and
    [PR #140](https://github.com/Parallel-NetCDF/PnetCDF/pull/140).

* Bug fixes
  +

* New example programs
  + C/flexible_bottom.c and C/vard_bottom.c - These two examples construct MPI
    derived data types using absolute memory addresses first and use MPI_BOTTOM
    when calling the PnetCDF flexible APIs.

* New programs for I/O benchmarks
  + C/pnetcdf_put_vara.c --
    * This program writes a series of 3D variables with 2D block-block
      partitioning pattern. Each variable is a record variable.
      See [PR #150](https://github.com/Parallel-NetCDF/PnetCDF/pull/150).
  + C/netcdf_put_vara.c --
    * This sequential NetCDF-C program writes a series of 3D variables. Each
      variable is a record variable.
    * This program and `C/pnetcdf_put_vara.c` can be used to compare the
      performance of NetCDF and PnetCDF when running sequentially, i.e. one
      process.
      See [PR #150](https://github.com/Parallel-NetCDF/PnetCDF/pull/150).

* New test program
  + test/testcases/flexible_large_count.c - tests flexible APIs that use MPI
    derived datatypes created by MPI large-count datatype constructors.
    See [PR #145](https://github.com/Parallel-NetCDF/PnetCDF/pull/145).

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
  + none

