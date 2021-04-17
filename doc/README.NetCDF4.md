# Support NetCDF-4 files

Starting from version 1.11.0, PnetCDF supports data access to HDF5-based
NetCDF-4 files. This format is also referred as the `enhanced data model`,
versus the `classic model` (in NetCDF language). The new I/O driver, called the
NetCDF-4 driver, is a wrapper of the NetCDF Library that translates PnetCDF
APIs to their corresponding NetCDF APIs, for example, `ncmpi_create` calls
`nc_create`. Through calling PnetCDF APIs, this feature allows applications to
read/write files in NetCDF-4 format.

## Enable NetCDF-4 support

PnetCDF requires a NetCDF library that is built with parallel HDF5 capability.
Example build instructions for HDF5 and NetCDF-4 are given below.
* To build parallel HDF5 (version 1.10.4 or later is recommended)
  + Source tar ball of HDF5 can be downloaded from URL:
    https://www.hdfgroup.org/downloads/hdf5/source-code/
  + Build commands:
    ```
    gzip -dc hdf5-1.10.4.tar.gz | tar -xf -
    cd hdf5-1.10.4
    ./configure --prefix=/HDF5/install/path \
                --enable-parallel=yes \
                CC=mpicc FC=mpifort CXX=mpicxx
    make install
    ```
* To build NetCDF-4 (version 4.6.2 or later is recommended)
  + Source tar ball of NetCDF-4 can be downloaded from URL:
    https://github.com/Unidata/netcdf-c/releases
  + Build commands:
    ```
    gzip -dc v4.6.2.tar.gz | tar -xf -
    cd netcdf-c-4.6.2
    ./configure --prefix=/NetCDF4/install/path \
                --enable-netcdf-4 \
                CC=mpicc \
                CPPFLAGS="-I/HDF5/install/path/include" \
                LDFLAGS="-L/HDF5/install/path/lib"
    make install
    ```
  + Note NetCDF-4 library built with PnetCDF feature enabled, i.e. option
    `--enable-pnetcf`, is not supported.
* To build PnetCDF with NetCDF-4 support
  + Add `--with-netcdf4` option at the configure command line. Option
    `--with-netcdf` can also be used to specify the installation path of
    NetCDF-4 if it is not in the standard localtions. For example,
    ```
    gzip -dc pnetcdf-1.11.0.tar.gz | tar -xf -
    cd pnetcdf-1.11.0
    ./configure --prefix=/PnetCDF/install/path \
                --with-netcdf4=/NetCDF4/install/path
    ```

## Accessing NetCDF-4 file

To create a NetCDF-4 file, add the flag NC_NETCDF4 into argument cmode when
calling `ncmpi_create()`. For example,
```
int cmode;
cmode = NC_CLOBBER | NC_NETCDF4;
ncmpi_create(MPI_COMM_WORLD, "testfile.nc", cmode, MPI_INFO_NULL, &ncid);
```
or
```
cmode = NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL;
ncmpi_create(MPI_COMM_WORLD, "testfile.nc", cmode, MPI_INFO_NULL, &ncid);
```

NC_NETCDF4 flag is not required when opening an existing file in NetCDF-4
format, as PnetCDF checks the file format and selects the proper I/O driver.

Users can also set the default file format to NetCDF-4 by calling API
`ncmpi_set_default_format()` with argument `NC_FORMAT_NETCDF4C` or
`NC_FORMAT_NETCDF4_CLASSIC`. For example,
```
ncmpi_set_default_format(NC_FORMAT_NETCDF4, &old_formatp);
```
or
```
ncmpi_set_default_format(NC_FORMAT_NETCDF4_CLASSIC, &old_formatp);
```
When no file format specific flag is set in argument cmode, PnetCDF will use
the default setting.


## Known Problems

Some features are not supported due to the availability of APIs different
between PnetCDF and NetCDF libraries. I/O semantics are also slightly
different.

* API families of `vard`, `varn`, and nonblocking I/O are not supported. This
  is because NetCDF does not have corresponding APIs. Error code
  `NC_ENOTSUPPORT` will be returned. For vard APIs, NetCDF does not have APIs
  that allow accessing the file directly by an MPI derived file type. For
  `varn` APIs, a potential solution is to split the request into into multiple
  vara calls. However, such solution must deal with the situation when the
  numbers of requests are different among processes in the collective data
  mode. Supporting `varn` is thus a future work.
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
  the HDF bug issue HDFFV-10501. The bug fix will appear in HDF5 1.10.4 release.
* PnetCDF allows different kinds of APIs called by different processes in a
  collective I/O. For example, `ncmpi_put_vara_int_all()` is called at process
  0 and `ncmpi_put_vars_float_int()` is called at process 1. However, this is
  not allowed in NetCDF-4. Thus, the same kind of API must be used in a
  collective call when accessing a NetCDF-4 file.
* When creating a new NetCDF-4 file, PnetCDF only allows creating 1 unlimited
  dimension. When reading an existing NetCDF-4 file, the file can have more
  than 1 unlimited dimension. However, in this case, API ncmpi_inq_unlimdim
  will only return the dimension ID of first unlimited dimension. This behavior
  conforms with NetCDF-4 library.
* The following APIs are not supported yet. An error code NC_ENOTSUPPORT will
  be returned, if called.
  * `ncmpi_inq_header_size`, `ncmpi_inq_header_extent`
  * `ncmpi_inq_striping`
  * `ncmpi_inq_varoffset`
  * `ncmpi_fill_var_rec`
  * `ncmpi_sync_numrecs`
  * `ncmpi_flush`
  * all nonblocking APIs
  * vard and varn APIs
  * flexible APIs (i.e. argument buftype is a constructed MPI derived data type)
* APIs `ncmpi_inq_get_size` and `ncmpi_inq_put_size` report the amount of data
  that has been read or written, but excluding the I/O to the file header. They
  are simply the size of data passed between PnetCDF and NetCDF, not the actual
  size read from or written to the file system.

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in top-level directory.

