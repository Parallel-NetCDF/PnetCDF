# Support NetCDF-4 files

Starting from version 1.11.0, PnetCDF supports data access to HDF5-based
NetCDF-4 files. This format is also referred the enhanced data model, versus to
the classic model (in NetCDF language). The new I/O driver, called the NetCDF-4
driver, is a wrapper of the NetCDF Library that translates PnetCDF APIs to
their NetCDF APIs. It allows applications to call PnetCDF APIs to read/write
files in NetCDF-4 format.

## Enable NetCDF-4 support

PnetCDF requires a NetCDF library that is built with parallel HDF5 capability.
An example build instruction for HDF5 and NetCDF-4 is given below.
* To build parallel HDF5 (version 1.10.2 or later is recommended)
  + Source tar ball of HDF5 can be downloaded from URL:
    https://www.hdfgroup.org/downloads/hdf5/source-code/
  + Build commands:
    ```
    ./configure --prefix=/HDF5/install/path \
                --enable-parallel=yes \
                CC=mpicc FC=mpifort CXX=mpicxx
    make install
    ```
* To build NetCDF-4 (version 4.6.1 or later is recommended)
  + Source tar ball of NetCDF-4 can be downloaded from URL:
    https://github.com/Unidata/netcdf-c/releases
  + Build commands:
    ```
    ./configure --prefix=/NetCDF4/install/path \
                --enable-netcdf-4 \
                CC=mpicc \
                CPPFLAGS="-I/HDF5/install/path/include" \
                LDFLAGS="-L/HDF5/install/path/lib"
    make install
    ```
* To build PnetCDF with NetCDF-4 support
  + Add `--enable-netcdf4` option at the configure command line. Option
    `--with-hdf5` can be used to specify the installation path of HDF5 and
    option `--with-netcdf` for the installation path of NetCDF-4. For example,
    ```
    ./configure --prefix=/PnetCDF/install/path \
                --enable-netcdf4 \
                --with-hdf5=/HDF5/install/path \
                --with-netcdf4=/NetCDF4/install/path
    ```

## Accessing NetCDF-4 file

To create a NetCDF-4 file, add the flag NC_NETCDF4 into argument cmode when
calling `ncmpi_create()`. For example,
```
int cmode;
cmode = NC_CLOBBER | NC_NETCDF4;
or
cmode = NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL;
ncmpi_create(MPI_COMM_WORLD, "testfile.nc", cmode, MPI_INFO_NULL>, &ncid);
```

NC_NETCDF4 flag is not required when opening an existing file in NetCDF-4
format, as PnetCDF checks the file format and selects the proper I/O driver.

Users can also set the default file format to NetCDF-4 by calling API
`ncmpi_set_default_format()` with argument `NC_FORMAT_NETCDF4C` or
`NC_FORMAT_NETCDF4_CLASSIC`. For example,
```
ncmpi_set_default_format(NC_FORMAT_NETCDF4, &old_formatp);
or
ncmpi_set_default_format(NC_FORMAT_NETCDF4_CLASSIC, &old_formatp);
```
When no file format specific flag is set in argument cmode, PnetCDF will use
the default setting.


## Known Problems

Some features are not supported due to the availability different of API between
PnetCDF and NetCDF libraries. I/O semantics are also slightly different.

* API families of `vard` and `varn` are not supported. This is because NetCDF
  does not have corresponding APIs. Error code `NC_ENOTSUPPORT` will be
  returned. For vard APIs, NetCDF does not have APIs that allow accessing the
  file directly by an MPI derived file type. For `varn` APIs, a potential
  solution is to split the request into into multiple vara calls. However, such
  solution must deal with the situation when the numbers of requests are
  different among processes in the collective data mode. Supporting `varn` is
  thus a future work.
* Although a new file of format NC_FORMAT_NETCDF4 or NC_FORMAT_NETCDF4_CLASSIC
  can be created, the I/O operations are limited to the NetCDF classic model
  I/O, This is because PnetCDF APIs do not include those for operating enhanced
  data objects, such as groups, compound data types, compression etc. As for
  reading a NetCDF-4 file created by NetCDF-4, PnetCDF supports the classic way
  of reading variables, attributes, and dimensions. Inquiry APIs for enhanced
  metadata is currently not supported.
* Due to a bug in HDF5 version 1.10.2 and prior, collective I/O on record
  variables with a subset of participating processes making zero-length
  requests may cause an HDF5 error and the program to hang. Readers refer to
  the HDF bug issue HDFFV-10501.
* PnetCDF allows different kinds of APIs called by different processes in a
  collective I/O. For example, `ncmpi_put_vara_int_all()` is called at process
  0 and `ncmpi_put_vars_float_int()` is called at process 1. However, this is
  not allowed in NetCDF-4. Thus, the same kind of API must be used in a
  collective call when accessing a NetCDF-4 file.

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in top-level directory.

