# Support variable chunking and compression

PnetCDF contains an experimental variable chunking and compression feature 
for classic NetCDF files.

For details about its design and implementation, please refer to:
Hou, Kaiyuan, et al. "Supporting Data Compression in PnetCDF." 
2021 IEEE International Conference on Big Data (Big Data). IEEE, 2021.

## Enable variable chunking support

* To build PnetCDF with variable chunking support
  + Add `--enable-chunking` option at the configure command line. For example,
    ```
    ./configure --prefix=/PnetCDF/install/path --enable-chunking
    ```
* To build deflate filter support for chunked variable
  + Add `--enable-zlib` option at the configure command line. Option
    `--with-zlib` can also be used to specify the installation path of
    zlib if it is not in the standard locations. For example,
    ```
    ./configure --prefix=/PnetCDF/install/path --enable-chunking --enable-zlib \
                --with-zlib=/zlib/install/path
    ```
* To build sz filter support for chunked variable
  + Add `--enable-sz` option at the configure command line. Option
    `--with-sz` can also be used to specify the installation path of
    sz if it is not in the standard locations. For example,
    ```
    ./configure --prefix=/PnetCDF/install/path --enable-chunking --enable-sz \
                --with-sz=/sz/install/path
    ```

## Enable variable chunking

To enable chunked storage layout for variables, set the file info "nc_chunking" 
to "enable". The chunking feature requires 64-bit NetCDF format (CDF5). 
For example,
```
    MPI_Info_create(&info);
    ncmpi_create(MPI_COMM_WORLD, fname, NC_64BIT_DATA, info, &ncid);
```
Alternatively, the file info can be set through the environment variable 
"PNETCDF_HINTS".
```
export PNETCDF_HINTS="nc_chunking=enable"
```
When chunking is enabled, all non-scalar variables will be stored in a chunked 
storage layout. Scalar variables are not chunked.

Users can also set the default filter for chunked variables. For example,
```
    MPI_Info_set(info, "nc_chunk_default_filter", "zlib");
```
or
```
export PNETCDF_HINTS="nc_chunking=enable;nc_chunk_default_filter=zlib"
```
The available filter options are none (default), zlib (deflate), sz.

## Define chunk dimension of variables

Applications can use the following APIs to set and get the chunk dimension of 
a variable.
```
    int ncmpi_var_set_chunk (int ncid, int varid, int *chunk_dim);
    int ncmpi_var_get_chunk (int ncid, int varid, int *chunk_dim);
```
For example:
```
    int dim[2] = {100, 100};
    int chunk_dim[2] = {10, 10};
    ncmpi_def_var (ncid, name, type, 2, dim, &varid)
    ncmpi_var_set_chunk (ncid, varid, chunk_dim);
```
For record variables, the chunk dimension along the record dimension is always 
1.
The default chunk dimension is the dimension of the variable except for the 
record dimension. By default, PnetCDF will create one chunk per record or 
variable.

## Define filter for chunked variables

Applications can use the following APIs to set and get the chunk dimension of 
a variable.
```
#define NC_FILTER_NONE  0
#define NC_FILTER_DEFLATE   2
#define NC_FILTER_SZ    3
int ncmpi_var_set_filter (int ncid, int varid, int filter);
int ncmpi_var_get_filter (int ncid, int varid, int *filter);
```
For example:
```
    ncmpi_var_set_filter (ncid, varid, NC_FILTER_DEFLATE);
```
Valid filter values are NC_FILTER_NONE (none), NC_FILTER_DEFLATE (zlib), and 
NC_FILTER_SZ (sz).


## Known problems

There are some limitations of the experimental variable chunking feature.

* Only one filter can be applied to a chunked variable. Unlike HDF5 which allows 
  the stacking of multiple filters on chunked datasets, the current 
  implementation in PnetCDF only allows a single filter to be applied to a 
  variable.
* No per-variable option for variable chunking. If chunking is enabled, all 
  non-scalar variables will be chunked even if the chunk dimension is not 
  defined.
* Independent variable I/O is not supported. Variable read/write (get/put) 
  must be collective in order to maintain data consistency of filtered chunks. 
  Non-blocking APIs can be used to mitigate the impact of this limitation.

Copyright (C) 2022, Northwestern University and Argonne National Laboratory

See the COPYRIGHT notice in the top-level directory.

