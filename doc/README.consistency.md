## Note on parallel I/O data consistency

PnetCDF follows the same parallel I/O data consistency as MPI-IO standard.
Refer the URL below for more information.
http://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node296.htm#Node296

Readers are also referred to the following paper.
Rajeev Thakur, William Gropp, and Ewing Lusk, On Implementing MPI-IO Portably
and with High Performance, in the Proceedings of the 6th Workshop on I/O in
Parallel and Distributed Systems, pp. 23-32, May 1999.

If users would like PnetCDF to enforce a stronger consistency, they should add
NC_SHARE flag when open/create the file. By doing so, PnetCDF adds
MPI_File_sync() after each MPI I/O calls.
  * For PnetCDF collective APIs, an MPI_Barrier() will also be called right
    after MPI_File_sync().
  * For independent APIs, there is no need for calling MPI_Barrier().

Users are warned that the I/O performance when using NC_SHARE flag could become
significantly slower than not using it.

If NC_SHARE is not set, then users are responsible for their desired data
consistency. To enforce a stronger consistency, users can explicitly call
ncmpi_sync(). In ncmpi_sync(), MPI_File_sync() and MPI_Barrier() are called.

### Note on header consistency in memory and file
In data mode, changes to file header can happen in the following scenarios.
  1. Renaming variables, dimensions, or attributes, when users make calls to
     ncmpi_rename_var(), ncmpi_rename_dim(), or ncmpi_rename_att().
  2. Overwriting an existing attribute with equal or smaller length, when users
     make calls to ncmpi_put_att() and its family APIs.
  3. New records are created from put APIs.

For cases 1 and 2, PnetCDF requires those APIs that change names or attributes
be called collectively. Internally, PnetCDF calls ncmpii_write_header() to
write (overwrite) the file header before these APIs return. Therefore, the file
header in memory is sync-ed and updated in file. Note starting from version
1.4.0, these APIs if called in data mode (independent or collective), they must
be called collectively by all processes that open the file.

For case 3, change to the number of records (numrecs) can happen in collective
or independent data mode by put APIs.
  * When in collective mode, PnetCDF sync numrecs in memory among all processes
    followed by letting root write the new value to the file. This applies to
    all collective blocking APIs and nonblocking ncmpi_wait_all().

  * When in independent mode, PnetCDF updates numrecs in local memory and sets
    the dirty bit NC_NDIRTY to the file flag locally. Sync in memory and update
    to file will wait until a collective call to:
    + ncmpi_end_indepi_data()
    + ncmpi_sync_numrecs()
    + ncmpi_sync()
    + ncmpi_abort()
    + ncmpi_redef()
    + ncmpi_close()

Some facts in PnetCDF implementation:
1. The only part of file header that can be out-of-sync is numrecs. The
   out-of-sync numrecs can only be caused by independent APIs.

2. In collective data mode, header (including numrecs) is always sync-ed in
   memory and updated in file before the collective API returns.

3. Renaming is allowed in data mode. When called in data mode, they must be
   called by all processes. The file header will be updated before the API
   returns.

4. ncmpi_put_att() and its family are allowed in data mode. When called in data
   mode, they must be called by all processes. The file header will be updated
   before the API returns.

5. rename() and put_att() are not allowed to be called independently in data
   mode.

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

