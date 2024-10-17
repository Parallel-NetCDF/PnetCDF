------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New features
  + Support MPI derived data types that are constructed from MPI 4.0 large-count
    derived datatype constructors.
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
  + none

* New constants
  + C macro NC_FillValue replaces _FillValue. This conforms with NetCDF4's
    change in its version 4.9.3 release.
    See [PR #153](https://github.com/Parallel-NetCDF/PnetCDF/pull/153).

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
  + Fix `pnetcdf-config` of reflecting the installation path when installation
    is done by running command `make install DESTDIR=/alternate/directory`
    which prepends '/alternate/directory' before all installation names.
    See [PR #154](https://github.com/Parallel-NetCDF/PnetCDF/pull/154).

* New example programs
  + C/flexible_bottom.c and C/vard_bottom.c - These two examples construct MPI
    derived data types using absolute memory addresses first and use MPI_BOTTOM
    when calling the PnetCDF flexible APIs.

* New programs for I/O benchmarks
  + C/pnetcdf_put_vara.c --
    * This program writes a series of 3D variables with 2D block-block
      partitioning pattern. Each variable is a record variable.

  + C/netcdf_put_vara.c --
    * This sequential NetCDF-C program writes a series of 3D variables. Each
      variable is a record variable.
    * This program and `C/pnetcdf_put_vara.c` can be used to compare the
      performance of NetCDF and PnetCDF when running sequentially, i.e. one
      process.

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

