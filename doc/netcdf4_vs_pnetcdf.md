# Comparison of PnetCDF and NetCDF4

* [Supported File Formats](#supported-file-formats)
* [Programming Differences](#programming-differences)
* [Define Mode and Data Mode](#define-mode-and-data-mode)
* [Collective and Independent I/O Mode](#collective-and-independent-io-mode)
* [Blocking vs. Nonblocking APIs](#blocking-vs-nonblocking-apis)

---

## Supported File Formats
* NetCDF4 supports both classic and HDF5-based file formats.
  + Classic file format (CDF-1) -- The ESDS Community Standard defined the file format
    to be used in the NetCDF user community in 1989. The file header bears a
    signature of character string 'CDF-1' and now is commonly referred to as
    [CDF-1](https://parallel-netcdf.github.io/doc/c-reference/pnetcdf-c/CDF_002d1-file-format-specification.html)
    file format.
    * 'CDF-2' format -- The CDF-1 format was later extended to support large
      file size (i.e.  larger than 2GB) in 2004. See its specification in
      [ESDS-RFC-011v2.0](https://cdn.earthdata.nasa.gov/conduit/upload/496/ESDS-RFC-011v2.00.pdf).
      Because its file header bears a signature of 'CDF-2' and the format is
      also commonly referred to as
      [CDF-2](https://parallel-netcdf.github.io/doc/c-reference/pnetcdf-c/CDF_002d2-file-format-specification.html)
      format.
    * [CDF-5](https://parallel-netcdf.github.io/doc/c-reference/pnetcdf-c/CDF_002d5-file-format-specification.html)
      format -- The CDF-2 format was extended by PnetCDF developer team
      in 2009 to support large variables and additional large data types, such
      as 64-bit integer.
  + HDF5-based file format -- Starting from its version 4.0.0, NetCDF includes
    the format that is based on HDF5, which is referred to as NetCDF-4 format.
    This offer new features such as groups, compound types, variable length
    arrays, new unsigned integer types, etc.
* PnetCDF supports only the classic file formats.
  + The classic files created by applications using NetCDF4 library can be read
    by the PnetCDF library and vice versa.
  + PnetCDF provides parallel I/O for accessing files in the classic format.
    NetCDF4's parallel I/O for classic files makes use of PnetCDF library
    underneath. Such feature can be enabled when building NetCDF4 library.


---

## Programming Differences
* The API names are different between NetCDF4 and PnetCDF.
  + For C programming, NetCDF4 uses prefix `nc_` while PnetCDF uses `ncmpi_`.
  + For Fortran 77 programming, NetCDF4 uses prefix `nf_` while PnetCDF uses `nfmpi_`.
  + For Fortran 90 programming, NetCDF4 uses prefix `nf90_` while PnetCDF uses `nf90mpi_`.
* The argument list of APIs are the same between NetCDF4 and PnetCDF, except for file open and create.
* Table below shows examples of source code.
  + Both example codes write to a variable named `WIND` in the collective mode.
  + The differences are marked in colors, green for NetCDF and blue for PnetCDF.

| NetCDF4 | PnetCDF |
|:-------|:--------|
| /* mode to create a file in classic CDF-5 format */ | |
| int cmode = NC_64BIT_DATA; | |
| /* create a new file */ | |
| ${\textsf{\color{green}nc\\_create\\_par}}$(path, cmode, comm, info, &ncid); | ${\textsf{\color{blue}ncmpi\\_create}}$(comm, path, cmode, info, &ncid);|
| /* add a global attributes */ | |
| char *attr = "Wed Mar 27 14:35:25 CDT 2024"; | |
| ${\textsf{\color{green}nc\\_put\\_att\\_text}}$(ncid, NC_GLOBAL, "history", 28, attr); | ${\textsf{\color{blue}ncmpi\\_put\\_att\\_text}}$(ncid, NC_GLOBAL, "history", 28, attr);  |
| /* define dimensions */ | |
| ${\textsf{\color{green}nc\\_def\\_dim}}$(ncid, "time", NC_UNLIMITED, &dimid[0]); | ${\textsf{\color{blue}ncmpi\\_def\\_dim}}$(ncid, "time", NC_UNLIMITED, &dimid[0]); |
| ${\textsf{\color{green}nc\\_def\\_dim}}$(ncid, "lat",  360,          &dimid[1]); | ${\textsf{\color{blue}ncmpi\\_def\\_dim}}$(ncid, "lat",  360,          &dimid[1]); |
| ${\textsf{\color{green}nc\\_def\\_dim}}$(ncid, "lon",  720,          &dimid[2]); | ${\textsf{\color{blue}ncmpi\\_def\\_dim}}$(ncid, "lon",  720,          &dimid[2]); |
| /* define a 3D variable of float type */ | |
| ${\textsf{\color{green}nc\\_def\\_var}}$(ncid, "WIND", NC_FLOAT, 3, dimid, &varid); | ${\textsf{\color{blue}ncmpi\\_def\\_var}}$(ncid, "WIND", NC_FLOAT, 3, dimid, &varid); |
| /* add attributes to the variable */ | |
| attr = "atmospheric wind velocity magnitude"; | |
| ${\textsf{\color{green}nc\\_put\\_att\\_text}}$(ncid, varid, "long_name", 35, attr); |${\textsf{\color{blue}ncmpi\\_put\\_att\\_text}}$(ncid, varid, "long_name", 35, attr); |
| /* exit define mode */ | |
| ${\textsf{\color{green}nc\\_enddef}}$(ncid); | ${\textsf{\color{blue}ncmpi\\_enddef}}$(ncid); | |
| /* collectively write to variable WIND */ | |
| ${\textsf{\color{green}size\\_t}}$ start[3], count[3]; | ${\textsf{\color{blue}MPI\\_Offset}}$ start[3], count[3]; |
| ${\textsf{\color{green}nc\\_var\\_par\\_access}}$(ncid, NC_GLOBAL, NC_COLLECTIVE); | |
| ${\textsf{\color{green}nc\\_put\\_vara\\_float}}$(ncid, varid, start, count, buf); | ${\textsf{\color{blue}ncmpi\\_put\\_vara\\_float\\_all}}$(ncid, varid, start, count, buf); |
| /* close file */ | |
| ${\textsf{\color{green}nc\\_close}}$(ncid); | ${\textsf{\color{blue}ncmpi\\_close}}$(ncid); |


---

## Define Mode and Data Mode

In PnetCDF, an opened file is in either define mode or data mode. Switching
between the modes is done by explicitly calling `"ncmpi_enddef()"` and
`"ncmpi_redef()"`. NetCDF4 when operating on an HDF5-based file has no such
mode switching requirement. The reason of PnetCDF enforcing such a requirement
is to ensure the metadata consistency across all the MPI processes and keep the
overhead of metadata synchronization small.

* Define mode
  + When calling `"ncmpi_create()"` to create a new file, the file is
    automatically put in the define mode.  While in the define mode, the user
    program can create new dimensions, new variables, and netCDF attributes.
    Modification of these data objects' metadata can only be done when the file
    is in the define mode.
  + When opening an existing file, the opened file is automatically put in the
    data mode. To add or modify the metadata, the user program must call
    `"ncmpi_redef()"`.

* Data mode
  + Once the creation or modification of metadata is complete, the user program
    must call `"ncmpi_enddef()"` to leave the define mode and enter the data
    mode.
  + While an open file is in data mode, the user program can make read and
    write requests to that variables that have been created.

<ul>
  <li> A PnetCDF example codes below show switching between define and data
       modes after creating a new file.</li>
  <li> <details>
  <summary>Example code fragment (click to expand)</summary>

```c
  #include <mpi.h>
  #include <pnetcdf.h>
  ...
  /* Create the file */
  ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);

  ...
  /* Define dimensions */
  ncmpi_def_dim(ncid, "Y", 16, &dimid[0]);
  ncmpi_def_dim(ncid, "X", 32, &dimid[1]);

  /* Define a 2D variable of integer type */
  ncmpi_def_var(ncid, "grid", NC_INT, 2, dimid, &varid);

  /* Add an attribute of string type to the variable */
  str_att = "example attribute of type text";
  ncmpi_put_att_text(ncid, varid, "str_att_name", strlen(str_att), str_att);

  /* Exit the define mode */
  ncmpi_enddef(ncid);

  /* Write to a subarray of the variable */
  MPI_Offset start[2], count[2];
  start[0] = 4;
  start[1] = 8;
  count[0] = 10;
  count[1] = 10;
  ncmpi_put_vara_int_all(ncid, varid, start, count, buf_int);

  /* Re-enter the define mode */
  ncmpi_redef(ncid);

  /* Define a new 2D variable of float type */
  ncmpi_def_var(ncid, "temperature", NC_FLOAT, 2, dimid, &var_flt);

  /* Exit the define mode */
  ncmpi_enddef(ncid);

  /* Write to a subarray of the variable, var_flt */
  start[0] = 2;
  start[1] = 8;
  count[0] = 5;
  count[1] = 5;
  ncmpi_put_vara_float_all(ncid, var_flt, start, count, buf_flt);

  /* Close the file */
  ncmpi_close(ncid);
```
</details></li>

  <li> An example shows switching between define and data modes after opening an existing file.
  </li>
  <li> <details>
  <summary>Example code fragment (click to expand)</summary>

```c
  #include <mpi.h>
  #include <pnetcdf.h>
  ...
  /* Opening an existing file */
  ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);

  ...
  /* get the ID of variable named 'grid', a 2D variable of integer type */
  ncmpi_inq_varid(ncid, "grid", &varid);

  /* Read the variable's attribute named "str_att_name" */
  char str_att[64];
  ncmpi_get_att_text(ncid, varid, "str_att_name", str_att);

  /* Read a subarray of the variable, var */
  MPI_Offset start[2], count[2];
  start[0] = 4;
  start[1] = 8;
  count[0] = 10;
  count[1] = 10;
  ncmpi_get_vara_int_all(ncid, varid, start, count, buf_int);

  /* Re-enter the define mode */
  ncmpi_redef(ncid);

  /* Define a new 2D variable of double type */
  ncmpi_def_var(ncid, "precipitation", NC_DOUBLE, 2, dimid, &var_dbl);

  /* Add an attribute of string type to the variable */
  str_att = "mm/s";
  ncmpi_put_att_text(ncid, var_dbl, "unit", strlen(str_att), str_att);

  /* Exit the define mode */
  ncmpi_enddef(ncid);

  /* Write to a subarray of the variable, var_dbl */
  start[0] = 2;
  start[1] = 8;
  count[0] = 5;
  count[1] = 5;
  ncmpi_put_vara_double_all(ncid, var_dbl, start, count, buf_dbl);

  /* Close the file */
  ncmpi_close(ncid);
```
</details></li>
</ul>


---
## Collective and Independent I/O Mode

The terminology of collective and independent I/O comes from MPI standard. A
collective I/O function call requires all the MPI processes opening the same
file to participate. On the other hand, an independent I/O function can be
called by an MPI process independently from others.

For metadata I/O, both PnetCDF and NetCDF4 require the function calls to be
collective.

* Mode Switch Mechanism
  + PnetCDF -- when a file is in the data mode, it can be put into either
    collective or independent I/O mode.  The default mode is collective I/O
    mode.  Switching to and exiting from the independent I/O mode is done by
    explicitly calling `"ncmpi_begin_indep_data()"` and
    `"ncmpi_end_indep_data()"`.

  + NetCDF4 -- collective and independent mode switching is done per variable
    basis. Switching mode is done by explicitly calling `"nc_var_par_access()"`
    before accessing the variable. For more information, see
    [Parallel I/O with NetCDF-4](https://docs.unidata.ucar.edu/netcdf-c/current/parallel_io.html).

<ul>
  <li> A PnetCDF example shows switching between collective and independent I/O
       modes.</li>
  <li> <details>
  <summary>Example code fragment (click to expand)</summary>

```c
  #include <mpi.h>
  #include <pnetcdf.h>
  ...
  /* Create the file */
  ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);

  ...
  /* Metadata operations to define dimensions and variables */
  ...
  /* Exit the define mode (by default, into the collective I/O mode) */
  ncmpi_enddef(ncid);

  /* Write to variables collectively */
  ncmpi_put_vara_int_all(ncid, varid, start, count, buf_int);

  /* Read from variables collectively */
  ncmpi_get_vara_float_all(ncid, var_flt, start, count, buf_flt);

  /* Leaving collective I/O mode and entering independent I/O mode */
  ncmpi_begin_indep_data(ncid);

  /* Write to variables independently */
  ncmpi_put_vara_int(ncid, varid, start, count, buf_int);

  /* Read from variables independently */
  ncmpi_get_vara_float(ncid, var_flt, start, count, buf_flt);

  /* Close the file */
  ncmpi_close(ncid);
```
</details></li>
</ul>

<ul>
  <li> A NetCDF4 example shows switching between collective and
       independent I/O modes.</li>
  <li> <details>
  <summary>Example code fragment (click to expand)</summary>

```c
  #include <mpi.h>
  #include <netcdf.h>
  #include <netcdf_par.h>
  ...
  /* Create the file */
  nc_create_par(filename, NC_CLOBBER, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);

  ...
  /* Metadata operations to define dimensions and variables */
  ...

  /* set the access method to use MPI collective I/O for all variables */
  nc_var_par_access(ncid, NC_GLOBAL, NC_COLLECTIVE);

  /* Write to variables collectively */
  nc_put_vara_int(ncid, varid, start, count, buf_int);

  /* Read from variables collectively */
  nc_get_vara_float(ncid, var_flt, start, count, buf_flt);

  /* set the access method to use MPI independent I/O for all variables */
  nc_var_par_access(ncid, NC_GLOBAL, NC_INDEPENDENT);

  /* Write to variables independently */
  nc_put_vara_int(ncid, varid, start, count, buf_int);

  /* Read from variables independently */
  nc_get_vara_float(ncid, var_flt, start, count, buf_flt);

  /* Close the file */
  nc_close(ncid);
```
</details></li>
</ul>

---

## Blocking vs Nonblocking APIs
* Blocking APIs -- All NetCDF4 APIs are blocking APIs. A blocking API means the
  call to the API will not return until the operation is completed. For
  example, a call to `nc_put_var_float()` will return only when the write data
  has been stored at the system space, e.g. file systems. Similarly, a call to
  `nc_get_var_float()` will only return when the user read buffer containing
  the data retrieved from the file. Therefore, when a series of `put/get`
  blocking APIs are called, these calls will be committed by the NetCDF4
  library one at a time, following the same order of the calls.
* Nonblocking APIs -- In addition to blocking APIs, PnetCDF provides the
  nonblocking version of the APIs. A nonblocking API means the call to the API
  will return as soon as the `put/get` request has been registered in the
  PnetCDF library. The commitment of the request may happen later, when a call
  to `ncmpi_wait_all/ncmpi_wait` is made. The nonblocking APIs are listed below.
  + `ncmpi_iput_var_xxx()` - posts a nonblocking request to write to a variable.
  + `ncmpi_iget_var_xxx()` - posts a nonblocking request to from from a variable.
  + `ncmpi_bput_var_xxx()` - posts a nonblocking, buffered request to write to a variable.
  + `ncmpi_iput_varn_xxx()` - posts a nonblocking request to write multiple subarrays to a variable.
  + `ncmpi_iget_varn_xxx()` - posts a nonblocking request to read multiple subarrays from a variable.
  + `ncmpi_bput_varn_xxx()` - posts a nonblocking, buffered request to write multiple subarrays to a variable.
  + `ncmpi_wait_all()` - waits for nonblocking requests to complete, using collective MPI-IO.
  + `ncmpi_wait()` - waits for nonblocking requests to complete, using independent MPI-IO.
  + `ncmpi_attach_buff()` - Let PnetCDF to allocate an internal buffer to cache bput write requests.
  + `File.detach_buff()` - Free the attached buffer.
* The advantage of using nonblocking APIs is when there are many small
  `put/get` requests and each of them has a small amount.  PnetCDF tries to
  aggregate and coalesce multiple registered nonblocking requests into a large
  one, because I/O usually performs better when the request amounts are large
  and contiguous. See example programs
  [nonblocking_write.c](../examples/C/nonblocking_write.c) and
  [bput_varn_int64.c](../examples/C/bput_varn_int64.c).
* Table below shows the difference in C programming between using blocking
  and nonblocking APIs.

| PnetCDF Blocking APIs | PnetCDF Nonblocking APIs |
|:-------|:--------|
| ...<br>/* define 3 variables of NC_FLOAT type */ ||
| ncmpi_def_var(ncid, "PSFC", NC_FLOAT, 2, dimid, &psfc);<br>ncmpi_def_var(ncid, "PRCP", NC_FLOAT, 2, dimid, &prcp);<br>ncmpi_def_var(ncid, "SNOW", NC_FLOAT, 2, dimid, &snow); | ditto |
| ... ||
| /* exit define mode and enter data mode */<br>ncmpi_enddef(ncid); | ditto |
| ...<br>/* Call blocking APIs to write 3 variables to the file */ | <br>/* Call nonblocking APIs to post 3 write requests */ |
| ncmpi_put_vara_float_all(ncid, psfc, start, count, buf_psfc);<br>ncmpi_put_vara_float_all(ncid, prcp, start, count, buf_prcp);<br> ncmpi_put_vara_float_all(ncid, snow, start, count, buf_snow);| ncmpi_iput_vara_float(ncid, psfc, start, count, buf_psfc, &req[0]);<br>ncmpi_iput_vara_float(ncid, prcp, start, count, buf_prcp, &req[1]);<br>ncmpi_iput_vara_float(ncid, snow, start, count, buf_snow, &req[2]);|
| | /* Wait for nonblocking requests to complete */<br>ncmpi_wait_all(3, reqs, errs)|


