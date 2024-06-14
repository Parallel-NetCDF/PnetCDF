## Programming difference between PnetCDF blocking APIs and nonblocking APIs

* **Blocking APIs** return only when the put/get operations is completed, i.e.
  the write data is safely stored in the file system, and the read data has
  been retrieved and stored in the user buffer.
* **Non-blocking APIs** return once the metadata of the put/get requests has
  been copied internally in PnetCDF. Each MPI process can post multiple,
  a different number of non-blocking put/get calls. All the pending
  non-blocking requests will be completed only when API
  `ncmpi_wait_all/ncmpi_wait` is called.
* Table below shows example source codes comparing usage of blocking and
  nonblocking APIs
  + The example codes write two variables named `WIND` and `RAIN`.
  + The code differences are marked in colors, green for blocking and blue for
    non-blocking.


| Blocking I/O | Non-blocking I/O |
|:-------|:--------|
| /* create a new file */ | |
| ncmpi_create(comm, path, mode, info, &ncid); | same |
| /* define dimensions */ | |
| ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dim[0]); | same |
| ncmpi_def_dim(ncid, "lat",  360,          &dim[1]); | same |
| ncmpi_def_dim(ncid, "lon",  720,          &dim[2]); | same |
| /* define a 3D variable of float type, named "WIND" */ | |
| ncmpi_def_var(ncid, "WIND", NC_FLOAT, 3, dim, &var[0]); | same |
| /* define a 3D variable of float type, named "RAIN" */ | |
| ncmpi_def_var(ncid, "RAIN", NC_FLOAT, 3, dim, &var[1]); | same |
| /* exit define mode */ | |
| ncmpi_enddef(ncid); | same |
| /* set subarray range accessed by this rank */ | |
| MPI_Offset start[3], count[3]; | same |
| start[0] = ...; start[1] = ...; start[2] = ...; | same |
| count[0] = ...; count[1] = ...; count[2] = ...; | same |
| | ${\textsf{\color{blue}int req[2], status[2];}}$ |
| ${\textsf{\color{green}/* collectively write to variable WIND */}}$ | ${\textsf{\color{blue}/* post a non-blocking write request */}}$ |
| ${\textsf{\color{green}ncmpi\\_put\\_vara\\_float\\_all}}$(ncid, var[0], start, count, windBuf); | ${\textsf{\color{blue}ncmpi\\_iput\\_vara\\_float}}$(ncid, var[0], start, count, windBuf, ${\textsf{\color{blue}\\&req[0]);}}$ |
| ${\textsf{\color{green}/* collectively write to variable RAIN */}}$ | ${\textsf{\color{blue}/* post a non-blocking write request */}}$ |
| ${\textsf{\color{green}ncmpi\\_put\\_vara\\_float\\_all}}$(ncid, var[1], start, count, rainBuf); | ${\textsf{\color{blue}ncmpi\\_iput\\_vara\\_float}}$(ncid, var[1], start, count, rainBuf, ${\textsf{\color{blue}\\&req[1]);}}$ |
| | ${\textsf{\color{blue}/* wait for non-blocking request to complete */}}$ |
| | ${\textsf{\color{blue}ncmpi\\_wait\\_all(ncid, 2, req, status);}}$ |
| /* close file */ | |
| ncmpi_close(ncid); | same |


