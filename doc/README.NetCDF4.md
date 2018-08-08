# Support NetCDF-4 files

PnetCDF supports data access to HDF5-based NetCDF-4 files. This format is also
referred the enhanced data model, versus to the classic model. The new I/O
driver, called the NetCDF-4 driver, is a wrapper of the NetCDF Library that
translates PnetCDF APIs to their NetCDF APIs. It allows to call PnetCDF APIs to
read/write files in NetCDF-4 format.

## Enable NetCDF-4 support

PnetCDF requires a NetCDF library that is built with parallel HDF5 capability.
An example build instruction for HDF5 and NetCDF-4 is given below.
* Build parallel HDF5
  + Source tar ball of HDF5 can be downloaded from URL:
    https://www.hdfgroup.org/downloads/hdf5/source-code/
  + Build commands:
    ```
    ./configure --prefix=/HDF5/install/path \
                --enable-parallel=yes \
                CC=mpicc FC=mpifort CXX=mpicxx
    make install
    ```
* Build NetCDF-4
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
int cmode = NC_CLOBBER | NC_NETCDF4;
ncmpi_create(MPI_COMM_WORLD, "testfile.nc", cmode, MPI_INFO_NULL>, &ncid);
```

NC_NETCDF4 flag is not required when opening an existing file in NetCDF-4
format, as PnetCDF checks the file format and selects the proper I/O driver.

Users can also set the default file format to NetCDF-4 by calling
`ncmpi_set_default_format()` using argument `NC_FORMAT_NETCDF4C`. For example,
```
ncmpi_set_default_format(NC_FORMAT_NETCDF4, &old_formatp);
```
When no file format specific flag is set in argument cmode, PnetCDF will use
the default.


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
* For writing, only the classic model `NC_FORMAT_NETCDF4_CLASSIC` is supported.
  This is mainly because the current PnetCDF APIs do not include those for
  operating enhanced data objects, such as groups, compound data types,
  compression etc. As for reading a NetCDF-4 file created by NetCDF-4, PnetCDF
  supports the classic way of reading variables, attributes, and dimensions.
  Inquiry APIs for enhanced metadata is currently not supported.
* Due to a bug in HDF5 (1.10.2 and prior at the time of this writing),
  collective I/O on variables with zero-length requests in some participating
  processes will cause an HDF5 error and the program may hang. Refer to the
  HDF issue HDFFV-10501.
* PnetCDF allows different kinds of APIs called by different processes in a
  collective I/O. For example, `ncmpi_put_vara_int_all()` is called at process
  0 and `ncmpi_put_vars_float_int()` is called at process 1. However, this is
  not allowed in NetCDF-4. Thus, the same kind of API must be used in a
  collective call when accessing a NetCDF-4 file.

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in top-level directory.

