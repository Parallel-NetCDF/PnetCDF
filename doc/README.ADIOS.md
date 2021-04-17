# Support Read Capability of BP files

Starting from version 1.12.0, PnetCDF supports read capability of BP files. This feature is implemented by adding an I/O driver that makes calls to the ADIOS library. Note that write capability is not yet supported.

## Enable ADIOS support

PnetCDF requires an [ADIOS library](https://www.olcf.ornl.gov/center-projects/adios/) version 1.x. We recommend version 1.13.1.

* To build ADIOS library
  + Obtain the source tarball of ADIOS from URL: https://www.olcf.ornl.gov/center-projects/adios/
  + Build command
    ```
    gzip -dc adios-1.13.1.tar.gz | tar -xf -
    cd adios-1.13.1
    export MPICC=<MPI C compiler command> 
    export MPICXX=<MPI C++ compiler command> 
    export MPIFC=<MPI Fortran compiler command> 
    export CC=<C compiler command> 
    export CXX=<C++ compiler command> 
    export FC=<Fortran compiler command> 
    ./configure --prefix=/ADIOS/install/path
    make install
    ```
* To build PnetCDF with ADIOS support
  + Add `--with-adios` option at the configure command line. Option
    `--with-adios` can also be used to specify the installation path of ADIOS
    if it is not in the standard locations. For example,
    ```
    gzip -dc pnetcdf-1.12.0.tar.gz | tar -xf -
    cd pnetcdf-1.12.0
    ./configure --prefix=/PnetCDF/install/path \
                --with-adios=/ADIOS/install/path
    ```

* For detailed ADIOS configuration options, please refer to [ADIOS User Manual]( https://users.nccs.gov/~pnorbert/ADIOS-UsersManual-1.13.0.pdf)

## Open BP files

PnetCDF checks the BP file format automatically. There is no need to add any flag to argument `omode` when calling `ncmpi_open()`. For example,
```
ncmpi_open(MPI_COMM_WORLD, "testfile.bp", NC_NOWRITE, MPI_INFO_NULL, &ncid);
```
Calling `ncmpi_open` with `NC_WRITE` flag set in argument `omode` will result in error code `NC_ENOTSUPPORT` returned.

## Example programs

Example programs are available in folder `./examples/adios`. Brief descriptions for all example programs can be found in `./examples/README`.

## Dumping BP files

the ncmpidump utility can be used to dump the content of a BP formatted file.

* To dump a BP file, use the utility as when dumping a NetCDF file.
  For example,
  ```
  $ src/utils/ncmpidump/ncmpidump test/adios/arrays.bp -h
  netcdf arrays {
  // file format: ADIOS BP Ver. 3
  dimensions:
          NX = 10 ;
          NY = 100 ;
  variables:
          double var_double_2Darray(NX, NY) ;
          int var_int_1Darray(NX) ;
  }
  ```

## Design of the ADIOS driver

The ADIOS driver is a wrapper of ADIOS function calls that enables users to read BP files using PnetCDF APIs. For some APIs, additional operations are required in order to translate the BP data structures into NetCDF semantics.

* File open
  + Detect BP file format
    + Because BP file format contains no file signature like NetCDF and HDF formats, checking the metadata is required.
    + BP files use a file footer (the last 28 bytes) to store the file offsets pointing to the location of metadata. 
      The footer contains three 64-bit unsigned integers representing the file offsets of metadata tables for process groups, variables, and attributes, respectively. The last 32-bit of the footer stores the BP version information.
    + The ADIOS driver calls the ADIOS file open API and if the call returns without an error, the file is considered a valid BP format.
  + Read metadata
    + Once a valid BP file is opened, ADIOS driver reads all the metadata to retrieve the information about the variables and dimensions. These metadata are cached in memory for fast access.
  + Parsing metadata
    + Dimension objects -- One of the differences between BP and NetCDF formats is the dimension object.
      In NetCDF, dimensions are objects shared by multiple variables. To find the size of a variable, an application first inquires dimensions used to define the variable and then inquires the size of each dimension object.
      In the BP format, there is no special object defined for dimensions. The dimension of a variable can be either a constant or a reference to a scalar variable.
    + Although the same scalar variable can be used as a dimension to define multiple variables, ADIOS library does not expose such information directly to the users, such as dimension IDs in NetCDF semantics. Instead, the dimensions of a variable are presented by their sizes when inquired.
    + Make use of `bp2ncd` utility -- ADIOS library distributions come with a utility program named `bp2ncd` that converts BP files to NetCDF files. We make use of its source codes to parse the variable and dimension information.
    + One-file-per-process cases -- `bp2ncd` does not support the one-file-per-process configuration. `bp2ncd` cannot retrieve the metadata of dimension and variables from all files. To resolve this, the PnetCDF ADIOS driver creates a virtualized dimension object for every dimension of each variable. It considers all dimensions of a variable are uniquely for that variable.

* Reading variables
  + Reference to variables
    + ADIOS library refers a variable by its name in the read APIs, unlike NetCDF by the variable’s ID. As a result, the PnetCDF ADIOS driver maps the variable ID to its name when calling the ADIOS read APIs.
  + Unsupported API kinds
    + Due to the difference in the API syntax between ADIOS and PnetCDF, not all PnetCDF APIs can be translated to ADIOS function calls. Below are the APIs that are currently not supported.
    + vard: ADIOS does not support the use of MPI derived datatypes.
    + varn: Although it is possible to translate each subarray request into an ADIOS read call, it does not perform well. 
    + `ncmpi_cancel`: canceling already posted non-blocking requests is not supported. This is due to no corresponding mechanism in ADIOS library to track the status of nonblocking requests.
  + Data provenance
    + The BP files can store the file access history. Variables can have more than one timestep and each timestep corresponds to update to the variables. For variables with more than one timestep, the PnetCDF’s ADIOS driver adds the timestep as an additional dimension (the most significant one) to those variables.
    
  + Data type conversion
    + ADIOS library puts users the responsibility of data type conversion, and therefore its read APIs return data in its origin type. PnetCDF, on the other hand, supports automatic type conversion. (In the meantime, PnetCDF also checks data type overflow during the conversion.) The PnetCDF ADIOS driver performs the conversion by first reading the variable into a temporary buffer, and then converting the data into the desired type.

* Reading attributes
  + ADIOS library has a very similar set of APIs for reading attributes as PnetCDF, which makes the mapping of PnetCDF APIs to ADIOS attribute APIs relatively easy.
  + Attributes of a variable
    + In BP, attributes that are represented by a pathname. For example, attribute `X` of variable `Y` is represented as `/Y/X`. When reading a variable attribute, the PnetCDF’s ADIOS driver prepends the variable’s name to the path and calls the ADIOS functions.
  + Global attributes
    + The pathname does not include any variable names.

* Inquiring dimensions
  + As described above, the dimension metadata is constructed at file open and cached in memory. Dimension inquiry APIs will be fulfilled by accessing the cached metadata.
  + Unlimited dimensions
    + The PnetCDF ADIOS driver considers all dimensions fixed-size. The inquiring APIs for record dimension will always return results as if no unlimited dimension is defined in the file. Because PnetCDF only supports read capability for the BP files, this limitation should not cause problems.


## Known issues

The following features are not supported because the availability of APIs is different between PnetCDF and ADIOS libraries. Some are because of the differences in I/O semantics.

* PnetCDF APIs that create or modify a BP file is not supported.
* The PnetCDF burst buffer feature is not supported for BP files.
* Flexible and non-blocking APIs are not supported. This
  is because ADIOS does not have corresponding APIs. Error code
  `NC_ENOTSUPPORT` will be returned.
* API `vard` family is not supported. ADIOS library does not have APIs that allow accessing the file directly by an MPI derived file type. 
* API `varn` family is not supported.
* Inquiring the record dimension will return as if no record dimension is defined.
* When reading the one-file-per-process BP files, the dimensions in individual files will be translated into virtualized dimensions that do not have the relationship to any scalar variables.
* ADIOS does not support collective I/O on reading. As a result, PnetCDF will always perform independent read even if collective API is called.
* Opening some of the BP files in big-endian on a little-endian machine may cause an error (ex. output file of ADIOS examples/C/attributes/attributes_write). This is caused by a bug in ADIOS library before (and including) version 1.13.1. When linking against ADIOS version before and including 1.13.1, PnetCDF will disable tests related to the bug. For details about the bug, please refer to https://github.com/ornladios/ADIOS/issues/195

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in the top-level directory.
