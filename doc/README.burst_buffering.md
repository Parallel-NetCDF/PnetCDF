## Using Burst Buffers in PnetCDF
Burst buffer driver implements a log-based I/O aggregation for write requests.
It is designed to work on a wide range of burst buffer architecture.

### Publication
For detailed description of the implementation of the burst buffer feature,
please refer to:

* Kai-Yuan Hou, Reda Al-Bahrani, Esteban Rangel, Ankit Agrawal, Robert Latham,
  Robert Ross, Alok Choudhary, and Wei-keng Liao. Integration of Burst Buffer
  in High-Level Parallel I/O Library for Exascale Computing Era. In the
  Workshop on Parallel Data Storage & Data Intensive Scalable Computing
  Systems, held in conjunction with the International Conference for High
  Performance Computing, Networking, Storage and Analysis, November 2018.
  http://cucis.ece.northwestern.edu/publications/pdf/HAR18.pdf

### Build PnetCDF with burst buffer feature
Add "--enable-burst-buffering" to your configure command line, e.g.
```console
    ./configure --prefix=/path/to/install --enable-burst-buffering
```

### Running applications to make use of burst buffers
The burst buffer feature is enabled by setting the PnetCDF I/O hint,
`nc_burst_buf`, in an MPI info object and passing it to file creation and
opening, for instance by adding the following line in the MPI program.
```c
    MPI_Info_set(info, "nc_burst_buf", "enable");
```
The hint can also be set through the environment variable `PNETCDF_HINTS` at
the run time.
```console
    export PNETCDF_HINTS="nc_burst_buf=enable"
```

### PnetCDF I/O hints for burst buffer controls

Below is a list of supported hints.
```
Hint key                        Values          Default  Description
---------                       ------          -------  -----------
nc_burst_buf                    enable/disable  disable  Enabling/disabling
                                                         the burst buffering.
nc_burst_buf_dirname            <Valid POSIX    ./       Directory where log
                                 Directory>              files will be stored.
                                                         This is the path burst
                                                         buffer is mounted.
nc_burst_buf_del_on_close       enable/disable  enable   Whether or not the log
                                                         files should be
                                                         deleted after the
                                                         NetCDF file is closed.
                                                         Disabling allows other
                                                         programs to access the
                                                         file.
nc_burst_buf_flush_buffer_size  <integer>       0        Amount of memory per
                                                         MPI process that allows
                                                         PnetCDF to allocate to
                                                         flush the logged data.
                                                         The unit is in bytes.
                                                         0 means unlimited.
                                                         User must guarantee
                                                         that it is larger than
                                                         any individual I/O
                                                         requests.
```

### Example job script using DataWarp on Cori @NERSC
```console
#!/bin/bash
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:10:00
#SBATCH -o output.txt
#DW jobdw capacity=1289GiB access_mode=private type=scratch pool=sm_pool
#
export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${DW_JOB_PRIVATE};nc_burst_buf_del_on_close=disable"
srun -n 1 ./a.out
```
Note in this example, hint `nc_burst_buf_dirname` is set to the DataWarp path
automatically by the job scheduler SLURM. See more information about the
DataWarp usage in:
http://www.nersc.gov/users/computational-systems/cori/burst-buffer

### Burst buffering design in PnetCDF

The burst buffer driver is a wrapper driver of the ncmpio (MPI-IO) driver. All
variable write APIs are intercepted and their requests are saved (cached) in
the burst buffer. File header I/O proceeds with the ncmpio driver as usual,
i.e. directly write to the header stored on the destination file system. When
flushing data stored in the burst buffer, the driver retrieves the write
requests and combines them into few large write requests, thus to achieve a
better I/O performance.  When the flushing buffer size is not big enough to
accommodate all data cached in burst buffer, flushing will be done in multiple
rounds.

The data stored in the burst buffer is flushed (log-replayed) when:
1. the NetCDF file is closed,
2. a read request is made to a variable written before,
3. `ncmpi_wait`/`ncmpi_wait_all` is called
4. `ncmpi_flush` is called, or
5. `ncmpi_close` is called.


### Known issues

1. Burst buffering delays file writes until log-replay time. If an error occurs
   to an individual request, it will be reported at the flushing time and only
   the first error encountered will be reported.

2. Partial flushing is not supported. Any flushing calls will flush the entire
   cached data to the destination file system.  Thus, cancelling nonblocking
   write requests may result in getting the error code `NC_EFLUSHED`, which
   means it is too late to cancel as the requests have been flushed.

3. Sequential consistency is not guaranteed. The burst buffer driver does not
   consider the order of write requests when flushing.  As a result, if the
   application write to the same file location multiple times without flushing,
   the resulting NetCDF file can contain either value regardless the order the
   write requests were made. Users must call `ncmpi_flush` in order to ensure the
   desired consistency.  For example, after the first write to a variable, a
   flush must be explicitly called before the second write to the same
   variable.

### Metadata file format
The file format uses the Backus Normal Form (BNF) grammar notation.
```c
metadata_file	= header [entry ...]	// two sections: header and entry
header		= magic format big_endian is_external num_ranks rank_id max_ndims num_entries entry_begin basename proc_name
magic		= "PnetCDF" VERSION	// case sensitive 8-byte string
format		= CDF_MAGIC | HDF5_MAGIC | BP_MAGIC
big_endian	= TRUE | FALSE		// Endianness of values stored in this
					// metadata file. TRUE means big Endian
is_external	= TRUE | FALSE		// whether data saved in the data file is
					// already in external representation
					// TRUE means type-convert/byte-swap is
					// already done
num_ranks	= INT64			// number of MPI processes
rank_id		= INT64			// MPI process rank ID
max_ndims	= INT64			// max ndims of all variables, this can
					// be calculated at enddef() and makes
					// fix-sized entries
num_entries	= INT64			// number of write requests
entry_begin	= INT64			// starting file offset of entry section
basename	= name			// base file name supplied by users
proc_name	= name			// MPI processor name, obtained from MPI_Get_processor_name()
					// can be used like constructing cb_nodes in ROMIO
name		= name_len namestring
name_len	= INT32
namestring	= CHAR [CHAR ...] padding
padding		= <0, 1, 2, or 3 bytes to next 4-byte boundary>	// padding uses null (\x00) bytes
CDF_MAGIC	= 'C' 'D' 'F' VERSION ZERO	// 8-byte string
HDF5_MAGIC	= "\211HDF\r\n\032\n"		// 8-byte string
BP_MAGIC	= 'B' 'P' VERSION padding ZERO	// 8-byte string
VERSION		= '0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9'
TRUE		= ONE
FALSE		= ZERO
CHAR		= <8-bit byte>
INT32		= <32-bit signed integer, native representation>
INT64		= <64-bit signed integer, native representation>


entry		= esize api_kind itype varid ndims data_off data_len start count stride
esize		= INT64				// byte size of an entry
api_kind	= VAR | VAR1 | VARA | VARS | INT32	// PnetCDF API kinds, positive INT represent varn with the number itself represent number of request
itype		= TEXT | SCHAR | UCHAR | SHORT |
		  USHORT | INT | UINT  | FLOAT |
		  DOUBLE | LONGLONG | ULONGLONG	// in-memory buffer data types
varid		= INT32				// variable ID
ndims		= INT32				// variable's number of dimensions
data_off	= INT64				// request's starting offset in data file
data_len	= INT64				// length of request occupied in data file
start		= [INT64 ...]			// number of elements == ndims * num_requests
count		= [INT64 ...]			// number of elements == ndims * num_requests or 0 (if counts is NULL)
stride		= [INT64 ...]			// number of elements == ndims or 0 (if not vars)

VAR		= NEGONE
VAR1		= NEGTWO
VARA		= NEGTHREE
VARS		= NEGFOUR

TEXT		= ONE
SCHAR		= TWO
UCHAR		= THREE
SHORT		= FOUR
USHORT		= FIVE
INT		= SIX
UINT		= SEVEN
FLOAT		= EIGHT
DOUBLE		= NINE
LONGLONG	= TEN
ULONGLONG	= ELEVEN
ZERO		= 0		// 4-byte integer in native representation
ONE		= 1		// 4-byte integer in native representation
TWO		= 2		// 4-byte integer in native representation
THREE		= 3		// 4-byte integer in native representation
FOUR		= 4		// 4-byte integer in native representation
FIVE		= 5		// 4-byte integer in native representation
SIX		= 6		// 4-byte integer in native representation
SEVEN		= 7		// 4-byte integer in native representation
EIGHT		= 8		// 4-byte integer in native representation
NINE		= 9		// 4-byte integer in native representation
TEN		= 10		// 4-byte integer in native representation
ELEVEN		= 11		// 4-byte integer in native representation
NEGONE		= -1		// 4-byte integer in native representation
NEGTWO		= -2		// 4-byte integer in native representation
NEGTHREE	= -3		// 4-byte integer in native representation
NEGFOUR		= -4		// 4-byte integer in native representation
```
* Metadata log file size:
  + header size = 80 bytes + strlen(basename) + strlen(proc_name)
  + each entry size = 48 bytes + ndims * 24 bytes

### data file format
```
data_file	= magic [entry ...]
magic		= "PnetCDF" VERSION	// case sensitive 8-byte string
entry		= <data>
```

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

