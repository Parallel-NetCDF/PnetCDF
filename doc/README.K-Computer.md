## Build PnetCDF on The K computer

This file contains information of building PnetCDF on the K computer at the
RIKEN Advanced Institute for Computational Science in Kobe, Japan.
http://www.aics.riken.jp/en/k-computer/about/
```
The MPI compilers on the K computer are wrappers of Fujitsu compilers:
    Language    MPI compiler    Sequential compiler
  -----------  --------------  --------------------
    Fortran     mpifrtpx        (frtpx)
    C           mpifccpx        (fccpx)
    C++         mpiFCCpx        (FCCpx)
```

The K computer is a cross-compile system.  Be sure to run configure with the
--build and --host flags to put it in "cross compile mode".  This will make
configure use compile-only tests, instead of the usual compile-and-run tests
(running tests on the K computer login node will not work.) Below shows a
configure command that was used to build PnetCDF successfully on K.
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
            MANIFEST_TOOL=: \
            TESTMPIRUN="mpiexec -n NP" \
            TESTSEQRUN="mpiexec -n 1"
```

For users who would customize the source codes and edit configure.ac, it is
necessary to run command "autoreconf" to generate an updated "configure" file.
However, there is a bug in autoconf 2.69 that causes building the PnetCDF
Fortran interfaces to fail. Please read file README.Fujitsu for the bug patch
to autoconf 2.69. If this bug fix is not applied to (or unable to due to
permission issue) the autoconf utility, please add LDFLAGS="-L." to the
configure command line.  Adding LDFLAGS="-L." is necessary as the Fujitsu
Fortran compiler requires no "-L" for its internal library files.  Without
this, configure command may fail.

To build shared library, add configure command-line option "--enable-shared"
and set the MPICC, MPICXX, and LDFLAGS environment variables to the followings.
```console
            MPICC="mpifccpx -Xg" \
            MPICXX="mpiFCCpx -Xg" \
            LDFLAGS="-Wl,-shared"
```

### Testing PnetCDF:

* One way to run "make check" or "make ptest" on K's compute nodes is through
  an interactive job.  An example command to submit an interactive job
  requesting 4 nodes for 20 minutes is shown below.
  ```console
    pjsub --interact --rsc-list "node=4" --rsc-list "elapse=00:20:00" \
          --sparam "wait-time=3600"
  ```
* Once in the interactive mode, run command
  ```console
    source /work/system/Env_base
    make -s check
    make -s ptest
  ```
* To exit the interactive mode, run command
  ```console
   exit
  ```

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

