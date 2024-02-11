## Note on parallel I/O data consistency

PnetCDF follows the same parallel I/O data consistency as MPI-IO standard,
quoted below.

```
Consistency semantics define the outcome of multiple accesses to a single file.
All file accesses in MPI are relative to a specific file handle created from a
collective open. MPI provides three levels of consistency:
  * sequential consistency among all accesses using a single file handle,
  * sequential consistency among all accesses using file handles created from a
    single collective open with atomic mode enabled, and
  * user-imposed consistency among accesses other than the above.
Sequential consistency means the behavior of a set of operations will be as if
the operations were performed in some serial order consistent with program
order; each access appears atomic, although the exact ordering of accesses is
unspecified. User-imposed consistency may be obtained using program order and
calls to MPI_FILE_SYNC.
```

Users are referred to the MPI standard Chapter 14.6 Consistency and Semantics
for more information.
http://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node296.htm#Node296

Readers are also referred to the following paper.
Rajeev Thakur, William Gropp, and Ewing Lusk, On Implementing MPI-IO Portably
and with High Performance, in the Proceedings of the 6th Workshop on I/O in
Parallel and Distributed Systems, pp. 23-32, May 1999.

* NC_SHARE has been deprecated in PnetCDF release of 1.13.0.
  + NC_SHARE is a legacy flag inherited from NetCDF-3, whose purpose is to
    provide some degree of data consistency for multiple processes concurrently
    accessing a shared file. To achieve a stronger consistency, user
    applications are required to also synchronize the processes, such as
    calling MPI_Barrier, together with nc_sync.
  + Because PnetCDF follows the MPI file consistency, which only addresses the
    case when all file accesses are relative to a specific file handle created
    from a collective open, NC_SHARE becomes invalid. Note that NetCDF-3
    supports only sequential I/O and thus has no collective file open per se.

If users would like a stronger consistency, they may consider using the code
fragment below after each collective/independent write API call (e.g.
`ncmpi_put_vara_int`, `ncmpi_put_vara_int_all`, `ncmpi_wait_all`
`ncmpi_enddef`, `ncmpi_redef`, `ncmpio_begin_indep_data`,
`ncmpio_end_indep_data`, etc.).
```
    ncmpi_sync(ncid);
    MPI_Barrier(comm);
    ncmpi_sync(ncid);
```
Users are warned that the I/O performance could become significantly slower.
Note `ncmpi_sync` is a collective call and can be called in either collective
or independent data mode.

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

