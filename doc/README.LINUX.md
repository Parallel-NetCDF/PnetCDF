## Linux cluster using Intel compilers

Building PnetCDF using Intel compilers are now requiring no special treatment.
Users are suggested to refer to file 'INSTALL' for instructions. Below is an
example configure command.
```console
./configure --prefix=/path/to/install \
            MPICC=mpiicc \
            MPICXX=mpiicpc \
            MPIF77=mpiifort \
            MPIF90=mpiifort
```

* Below are some instructions for computer platforms from earlier development
  time of PnetCDF. They can be safely ignored if you are building PnetCDF
  1.10.0 and later.

* here are the steps John Tannahill <tannahill1@llnl.gov> used when building
  PnetCDF on the "mcr" cluster at LLNL.  It is a Linux cluster with the Intel
  compilers.
```console
   setenv MPICC  mpiicc
   setenv MPIF77 mpiifc
   setenv F77    ifc
   setenv FC     ifc
   setenv CC     icc
```
* then run the usual "configure; make ; make install"

On 5 October 2005, Richard Hedges reported that FFLAGS and CFLAGS needed the -O
and -mp options or the test suite would report errors in the
n[c,f]mpi_put_var*_float routines:
```console
FFLAGS = -O -mp
CFLAGS = -O -mp
```

## Jazz, a Linux cluster @ANL using Intel compilers

```console
# ;cat ~/.soft
# #
# # This is the .soft file.
# # It is used to customize your environment by setting up environment
# # variables such as PATH and MANPATH.
# # To learn what can be in this file, use 'man softenv'.
# #
# #
# @default
# +pbs
# +intel-7.0
#
# ---
#
# ;which mpicc
# /soft/apps/packages/mpich-gm-1.2.5..9-pre6-gm-1.6.3-intel-7.0/bin/mpicc
# ;which mpif77
# /soft/apps/packages/mpich-gm-1.2.5..9-pre6-gm-1.6.3-intel-7.0/bin/mpif77
# ;which icc
# /soft/com/packages/intel-7/compiler70/ia32/bin/icc
# ;which ifc
# /soft/com/packages/intel-7/compiler70/ia32/bin/ifc
#
# ---

setenv CC icc
setenv FC ifc
setenv F90 ifc
setenv CXX icc
setenv MPICC mpicc
setenv MPIF77 mpif77
```

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

