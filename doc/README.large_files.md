## Note on Large-file Support

For latest information of large-file support, please visit PnetCDF home page:
https://parallel-netcdf.github.io and refer to the Section "A Note About Large
File Support".


### BACKGROUND
The classic NetCDF file format (CDF-1) uses a 4-byte value to hold the
offset into the file where one can find a variable's data

In October 2003, Greg Sjaardema <gdsjaar@sandia.gov> proposed a new file format
(CDF-2) which uses an 8-byte value for the offset. We use his approach in
PnetCDF, though we have modified his patch against netcdf-3.5-beta1 to apply to
our codebase.

I couldn't find a URL to Greg's original message, but here's Russ Rew's
followup:
http://www.unidata.ucar.edu/projects/coohl/mhonarc/MailArchives/netcdf/msg04811.html

This means there are two different but compatible implementations of the CDF-2
file format.  We (the PnetCDF developers) will make our best effort to keep our
implementation compatible with the serial netcdf implementation. Please report
any incompatibilities to the developers at parallel-netcdf@mcs.anl.gov.

### PRELIMINARIES
First, it is important that your MPI-IO implementation uses an 8 byte type to
represent the 'MPI_Offset'. Starting from 1.9.0, PnetCDF requires the
underneath MPI library's MPI_Offset to be of 8-byte size.

Second, your platform should use an 8 byte type to represent the 'off_t' type.
On Linux, solaris, IRIX64 (and quite possibly others), PnetCDF will
automatically add the right options to the compiler to make this happen.

Run configure as you normally would.  Let the developers know if configure says
your 'off_t' is 4 bytes.  Proceed to compile and install the library.

### USAGE
By default, PnetCDF will create CDF-1 formatted files. This will ensure that
datasets created by our library will be compatible with the large body of
applications which expect NetCDF files to be CDF-1 formatted.

To write a CDF-2 formatted file, add the flag 'NC_64BIT_OFFSET' to the
ncmpi_create() function call ( or nfmpi_create() if using the Fortran
interface)

The PnetCDF library will detect the format of the dataset. There are no
special options needed to read back files created with the
NC_64BIT_OFFSET or NC_64BIT_DATA flag set.

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

