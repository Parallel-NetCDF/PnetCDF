## Programming difference between NetCDF4 and PnetCDF

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
