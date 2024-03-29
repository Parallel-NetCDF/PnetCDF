## Note on Building PnetCDF using Fujitsu Compilers

This file is for users who will edit file "configure.in" and run command
"autoreconf" to regenerate file "configure".

Autoconf 2.69 has a problem of checking the version of the Fujitsu Fortran
compiler, which causes errors for building PnetCDF at the configure time.

The fix and patch have been discussed and provided in the link below. PnetCDF
users/developers are encouraged to apply to patch to their autoconf utility.

Subject: Support Fujitsu in _AC_PROG_FC_V
https://lists.gnu.org/archive/html/autoconf-patches/2013-11/msg00001.html

Once the patch is applied and autoconf rebuilt, PnetCDF can be built with
one of the following configure commands.

* for non-cross-compiling
  ```console
  ./configure --prefix=/path/to/install \
              --with-mpi=/path/to/MPI/installation
  ```
* for cross-compiling
  ```console
  ./configure --prefix=/path/to/install \
              --with-mpi=/path/to/MPI/installation
              cross_compiling="yes"
  ```
* for cross-compiling and wth explicit MPI compiler names
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=/path/to/MPI/installation/mpifccpx \
              MPICXX=/path/to/MPI/installation/mpiFCCpx \
              MPIF77=/path/to/MPI/installation/mpifrtpx \
              MPIF90=/path/to/MPI/installation/mpifrtpx
  ```
* for cross-compiling and wth explicit MPI compiler names and optimization
  flags
  ```console
  ./configure --prefix=/path/to/install \
              --host=sparc64-unknown-linux-gnu \
              --build=x86_64-redhat-linux \
              MPICC=mpifccpx \
              MPICXX=mpiFCCpx \
              MPIF77=mpifrtpx \
              MPIF90=mpifrtpx \
              CFLAGS="-Kfast" \
              CXXFLAGS="-Kfast" \
              FFLAGS="-Kfast" \
              FCFLAGS="-Kfast" \
              TESTMPIRUN="mpiexec -n NP" \
              TESTSEQRUN="mpiexec -n 1"
  ```


If the above patch is not (or unable to) applied to the autoconf utility,
please add LDFLAGS="-L." to the command line when running configure. This is
necessary to prevent the failure of building PnetCDF Fortran interfaces, due to
Fujitsu Fortran compiler requires no "-L" for linking internal libraries.
Adding "-L." to force the linker to use this flag.

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

