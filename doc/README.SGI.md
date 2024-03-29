### SGI UV 2000
* Endeavour @ NASA http://www.nas.nasa.gov/hecc/resources/endeavour.html
* There are 2 compilers available on Endeavour: Intel and GNU.  Intel compiler
  is recommended.
* To use Intel compiler, run command below to load the Intel compiler module.
  ```console
  module load comp-intel
  module load mpi-intel
  ```

* run configure command:
  ```console
  ./configure --prefix=/path/to/install \
              MPICC=icc MPICXX=icpc MPIF77=ifort MPIF90=ifort \
              CFLAGS="-O2" FFLAGS="-O2" FCFLAGS="-O2" \
              LIBS=-lmpi LDFLAGS=
  ```

* Some users have their environment variable LDFLAGS set to "-shared",
  which can prevent PnetCDF to build correctly. Please note PnetCDF
  currently supports to build a static library only.

* It is also possible to build PnetCDF with GNU-based MPI compilers on
  Endeavour. Here is configure command.
  ```console
  ./configure --prefix=/path/to/install \
              CFLAGS="-O2" FFLAGS="-O2" FCFLAGS="-O2" \
              LDFLAGS=
  ```

### SGI IRIX64
The PnetCDF library should build fine under IRIX.  There are just a few
issues not directly related to the library:

* The IRIX compiler does not like the cnap.c test.  gcc compiles the
  file without any warnings, even when the tests are configured with
  "--enable-strict".  Renier Vogelsang <reiner@sgi.com> reports that the
  latest MIPSpro compiler release fixes this issue, so upgrade if
  possible.

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

