# Support Read Capability of ADIOS BP files

Starting from version 1.12.0, PnetCDF supports read capability of BP files. Note that write capability is not yet supported in PnetCDF.

## Enable ADIOS support

PnetCDF requires an [ADIOS library](https://www.olcf.ornl.gov/center-projects/adios/) version 1.x. We recommend version 1.13.1.

* To build ADIOS library
  + Obtain the source tar ball of ADIOS from URL: https://www.olcf.ornl.gov/center-projects/adios/
  + Build command
    ```
    gzip -dc adios-1.13.1.tar.gz | tar –xf -
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

* For detailed ADIOS configuration options, please refer to [ADIOS User Manual]( https://users.nccs.gov/~pnorbert/ADIOS-UsersManual-1.13.0.pdf)

## Reading BP files

PnetCDF checks the BP file format automatically. There is no need to add any flag to argument `omode` when calling `ncmpi_open()`. For example,
```
ncmpi_open(MPI_COMM_WORLD, "testfile.bp", NC_NOWRITE, MPI_INFO_NULL, &ncid);
```
Calling `ncmpi_open` with `NC_WRITE` flag set in argument `omode` will result in error code `NC_ENOTSUPPORT` returned.

## Example programs

Example programs are available in folder `./examples/adios`. Brief descriptions for all example programs can be found in ` examples/README`.


## Design of ADIOS driver

The ADIOS driver is a wrapper of the ADIOS library that enables users to read BP files using PnetCDF APIs. The ADIOS driver maps ADIOS data structures to NetCDF data structures and translates PnetCDF API calls into corresponding ADIOS function calls.

* Opening the file
  When the application opens a BP file using PnetCDF, the library detects the BP file format and load the ADIOS driver to handle the file.
  + Determine file format
    + To our best knowledge, BP file format has no (specified) signature that can be used to recognize BP formatted files.
    + Instead of a header, the BP file format uses a footer to index the data. 
      The last 28 byte of the BP file is a mini footer that contains 3 64 bits unsigned integer that represents the position of process groups, variables, and attributes index table in the footer respectively followed by a 4 bytes version information. The only rule regarding the mini footer is that process group index comes before the variable index and variable index comes before attributes index.
    To determine the file format, the ADIOS driver first checks if the mini footer matches the BP specification. If so, the ADIOS driver tries to open the file using the ADIOS read API. If ADIOS file open does not report an error, it is considered a valid BP formatted file.
  When a BP file is opened by the ADIOS driver, it opens the file using ADIOS read API and keep the opened file handle. At the same time, the ADIOS driver parses the file metadata regarding variables and dimensions. Theses metadata is cached inside the driver to support efficient query from the application. 
  + Parsing the metadata
    + Lack of information regarding dimension
      A major difference between BP and NetCDF file is the handling of dimensions. In a NetCDF file, dimensions are a type of object that can be associated with multiple variables in the same file. To find out the size of a variable, the application query for dimensions associated with that variable and then query the size of each dimension object.
      On the other hand, there is no dimension object in the BP files. When defining variables, the application specifies its size along each dimension directly. The size can be specified either by a constant or a reference to a scalar variable where the value indicates the size. In the latter case, the variable being referred to serves a similar role of dimension object in the NetCDF file.
      Although BP file does allow linking variables' shape to scalar variables, ADIOS does not expose such relationship to the user. ADIOS always present the shape of a variable purely by its value even if it is defined by scalar variables. As a result, it is not possible to derive dimension information by calling ADIOS APIs. In fact, such information may not even exist in a BP file since variable shape can be defined purely by value without linking to any scalar variable.
    + Following the view of bp2ncd conversion utility
      ADIOS distribution comes with a NetCDF conversion utility that converts BP formatted file into NetCDF file called bp2ncd. We rely on the code in this utility to provide variable and dimension information. By using it, we not only get the dimension information we need but also ensures that the view of a BP file represented by PnetCDF is consistent with the NetCDF file converted by the utility program. We ran a modified version of bp2ncd in which it does not actually generate a NetCDF file, instead, it records every dimension and variable it tries to create. The ADIOS driver uses this record to infer dimension and variable information.     
    + Alternative approach
      The bp2ncd utility was not designed to convert all BP files. For example, the BP file that is stored in file-per-process fashion (POSIX transmission method) is not supported by bp2ncd. On those files, the utility code cannot be used to acquire metadata about dimension and variables. In such case, the ADIOS driver simply creates a virtualized dimension for every dimension of every variable. It assumes that all variable shapes are defined by constant value directly and no scalar variable is used as dimension size.
* Reading the variable
  The view of ADIOS driver on variables is based on the record of bp2ncd which differs to the view given by ADIOS read API. As a result, the ADIOS driver specifies variables by name instead of ID when reading variable data using the ADIOS read API.
  + Unsupported API type
    Due to the difference in the interface of ADIOS and PnetCDF, not all PnetCDF function call can be translated to ADIOS function call. Some of the functions are possible to translate to ADIOS but can not be done efficiently. Those functions are currently not supported by the ADIOS driver.
    + vard: ADIOS does not allow reading variable using derived datatype
    + low-level API: ADIOS does not allow reading variable using derived datatype
    + write related API: ADIOS driver is read-only
    + varn: Although it is possible to translate it into multiple ADIOS read call, it does not provide the performance advantage varn API is supposed to provide. 
    + non-blocking: Although ADIOS have a set of non-blocking APIs, it serves a different purpose to the non-blocking API in PnetCDF and hence cannot be used to support non-blocking API in PnetCDF.
  + Timestep
    The BP file not only record the current state of the dataset, but it also records the entire history. It introduces the concept of timestep to deal with different versions of the file. 
    Each timestep corresponds to a session between file opening and closing. Variables can have different values on different time steps.
    Since NetCDF file, as well as PnetCDF API, do not include such concept, the ADIOS driver always returns the last (latest) time step of a variable.
  + Type conversion
    ADIOS read API always return data in origin type associated with a variable or attribute. 
    PnetCDF, on the other hand, does allow users to specify a different datatype when reading a variable or attribute. 
    If the user specifies a different type than the variable or attribute, the ADIOS driver will read the variable into a temporary buffer, and then convert the data into the desired type by performing C type casting on each individual value.
* Reading Attributes
  Unlike variable and dimensions, metadata regarding dimensions are not parsed using the bp2ncd code at file opening stage. This is because ADIOS has a very similar interface regarding attributes to PnetCDF. PnetCDF attribute read request can be easily mapped to ADIOS attribute read request.
  + Unlike in PnetCDF where the user can choose to query attributes on a specific variable or on the file (global), ADIOS read API presents all attributes as global attributes (on the file). Attributes that are intended to be associated with a variable are indicated by giving it a path-like name with variable name as a prefix. For example, an attribute X associated with variable Y will be presented as a global named /Y/X.
  + When the user reads a global attribute, we pass it directly to ADIOS API. When the user reads a variable attribute, we append the name of the variable in question before the attribute name and pass to ADIOS read API. In this way, variable attributes can be accessed either globally using its full path or locally at the variable with the only attribute name.
* Querying Dimensions
  As mentioned in the open section, ADIOS API does not expose dimension information to the user. BP file format does not enforce explicit dimension declaration as well. 
  We rely on the cached metadata to respond to dimension related queries.
  + Unlimited Dimensions
    For performance consideration, we do not record variable write from the bp2ncd utility. Hence, the size of the unlimited dimension is not in the cached metadata. ADIOS read API also does not expose the size of the unlimited dimension. The ADIOS driver assumes the size of the unlimited dimension is equal to the largest it ever appears in variable sizes returned by ADIOS read API.  

## Known Problems

Some features are not supported because the availability of APIs is different
between PnetCDF and ADIOS libraries. I/O semantics are also slightly
different.

* The ADIOS driver is ready only. It does not support any API that creates or modifies 
  the BP file.
* The ADIOS driver does not work with PnetCDF burst buffer feature due to its ready only nature.
* API families of `vard`, `varn`, derived type, and nonblocking I/O are not supported. This
  is because ADIOS does not have corresponding APIs. Error code
  `NC_ENOTSUPPORT` will be returned. For vard APIs, ADIOS does not have APIs that allow accessing the file directly by an MPI derived file type. For
  `varn` APIs, a potential solution is to split the request into multiple vara calls. However, doing this will not provide any performance improvement over fragmented request as varn API is supposed to do.
* PnetCDF current does not recognize record dimension. Variable with record dimension can still be read, but PnetCDF will not return information regarding the number of records.
* Subfiled BP file is not fully supported. As mentioned in the design section, we cannot parse the dimension information in sub-filed BP file. In such case, the application will see a virtualized dimension that enables them to read variable data but does not provide the relationship of a variable to scalar variables (if any) that defines its shape.
* C memory layout (row major) is always assumed even calling from PnetCDF Fortran API.

Copyright (C) 2018, Northwestern University and Argonne National Laboratory

See COPYRIGHT notice in the top-level directory.

