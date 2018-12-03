# Support ADIOS BP files

Starting from version 1.12.0, PnetCDF supports reading data from ADIOS 1.x BP files. 
Through calling PnetCDF APIs, this feature allows applications to
read files in BP format. For now, PnetCDF cannot write to BP formated files.

## Enable ADIOS support

PnetCDF requires a ADIOS library that is built with parallel I/O support.
* To build PnetCDF with ADIOS support
  + Add `--enable-adios` option at the configure command line. Option
    `--with-adios` can be used to specify the installation path of ADIOS.
    For example,
    ```
    gzip -dc parallel-netcdf-1.11.0.tar.gz | tar -xf -
    cd parallel-netcdf-1.11.0
    ./configure --prefix=/PnetCDF/install/path \
                --enable-adios \
                --with-adios=/ADIOS/install/path
    ```

## Accessing ADIOS file

For now, PnetCDF can only read BP formated file. To open a AIODS BP file, add the flag NC_NOWRITE into argument cmode when
calling `ncmpi_open()`. For example,
```
int cmode;
cmode = NC_NOWRITE;
ncmpi_open(MPI_COMM_WORLD, "testfile.bp", cmode, MPI_INFO_NULL, &ncid);
```

Setting NC_WRITE flag will result in error. PnetCDF will recognize ADIOS BP file aumatically and selects the proper I/O driver.
No flag regarding file format is required.

## Example programs

An example program is avaiable at /examples/adios/read_all.c

## Known Problems

Some features are not supported due to the availability of APIs different
between PnetCDF and ADIOS libraries. I/O semantics are also slightly
different.

* The ADIOS driver is ready only driver. It does not support any API that creates or modifies 
  objects within the BP file.
* The ADIOS driver does not work with burst buffer driver due to its ready only nature.
* API families of `vard`, `varn`, derived type, and nonblocking I/O are not supported. This
  is because ADIOS does not have corresponding APIs. Error code
  `NC_ENOTSUPPORT` will be returned. For vard APIs, ADIOS does not have APIs
  that allow accessing the file directly by an MPI derived file type. For
  `varn` APIs, a potential solution is to split the request into into multiple
  vara calls. However, such solution must deal with the situation when the
  numbers of requests are different among processes in the collective data
  mode. Supporting `varn` is thus a future work.
* PnetCDF current does not recognize record dimension. Variable with record dimension can 
  still be read, but PnetCDF will not return information regarding number of records.

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in top-level directory.

