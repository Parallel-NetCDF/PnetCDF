## Build instructions for Cray XC40

* [Cori @ NERSC](http://www.nersc.gov/systems/cori/)
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn
  ```
  + Note that the default setting of PnetCDF is to build static libraries only.
    To also build shared libraries, add "--enable-shared" to the configure
    command.
    ```console
    ./configure --prefix=/path/to/install \
                cross_compiling="yes" \
                MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
                --enable-shared
    ```
  + Adding "-dynamic" to the environment variable LDFLAGS is no longer
    required when building shared libraries, as the default linking mode has
    been changed since CDT/19.06 in Nov. 2019. See Cori Programming Environment
    Change at
    https://docs.nersc.gov/systems/cori/timeline/default_PE_history/2019Dec-2020Jan/

* [Theta @ ALCF](https://www.alcf.anl.gov/theta)
  ```console
    ./configure --prefix=/path/to/install \
                cross_compiling="yes" \
                MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn
  ```
  + Note that the default setting of PnetCDF is to build static libraries only.
    To also build shared libraries, add "--enable-shared LDFLAGS=-dynamic" to
    the configure command.
    ```console
    ./configure --prefix=/path/to/install \
                cross_compiling="yes" \
                MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
                --enable-shared LDFLAGS=-dynamic
    ```
  + See Compiling and Linking Overview on Theta at
    https://www.alcf.anl.gov/support-center/theta/compiling-and-linking-overview-theta

* Note on running "make check" or "make ptest".
  Cori and Theta are in a cross-compile environment. Jobs must be run through a
  batch scheduler. At the time of this writing, the schedule is Slurm. Below
  configure command line adds 3 variables in preparation for job submission.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              TESTMPIRUN="srun -n NP" TESTSEQRUN="srun -n 1" \
              TESTOUTDIR="$SCRATCH"
  ```
  + Run test programs on a single MPI processes with command "make check".
    Add below lines to your batch script which runs a single MPI process on one
    node only. Variable TESTOUTDIR should be set to the file directory for
    storing the test output files.
    ```console
    #SBATCH -N 1
    #SBATCH -t 01:30:00
    make check
    ```
  + Run test programs in parallel. A sample batch script fragment shown below
    allocates 4 compute nodes. "make ptest" requires at least 4 MPI processes.
    The number can be changed to use less number of nodes. Variable TESTOUTDIR
    should be set to the file directory for storing the test output files.
    ```console
    #SBATCH -N 4
    #SBATCH -t 00:20:00
    make ptest
    ```
* Note on using GNU compilers (i.e. gcc, g++, and gfortran)
  To swap to the GNU programming environment from default, run command below
  first before running configure command:
  ```console
  module swap PrgEnv-intel PrgEnv-gnu
  ```

* Note on using Cray compilers (i.e. craycc, crayCC, crayftn)
  (Not supported on Theta. See note below.)
  To swap to the Cray programming environment from default, run command below
  first before running configure command:
  ```console
  module swap PrgEnv-intel PrgEnv-cray
  ```

### Historical notes before upgrade to CDT/19.06

* Building PnetCDF on the Cray XC40, tested on
  [Cori @ NERSC](http://www.nersc.gov/systems/cori) and
  [Theta @ ALCF](https://www.alcf.anl.gov/theta)
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              CFLAGS="" CXXFLAGS="" FFLAGS="" FCFLAGS=""
  ```

  On Cori, when using Intel-based compilers, NERSC recommends the default level
  of optimization, i.e. no optimization arguments to the compiler. See section
  "Recommended flags" for Intel compilers at
  http://www.nersc.gov/users/computational-systems/cori/programming/compiling-codes-on-cori/

  * On April 29, 2019, The default gcc module version 7.0.3 is loaded on Cori,
    which causes `make` process to fail with an error message of
    ```console
    /usr/bin/ld: attempted static link of dynamic object `/opt/gcc/7.3.0/snos/lib/../lib64/libstdc++.so'
    ```
    This problem can be resolved by unloading the new gcc module before running
    the configure command, i.e.
    ```console
    module unload gcc/7.3.0
    ```

* Note that the default setting of PnetCDF is to build static libraries only.
  To also build shared libraries, add "--enable-shared LDFLAGS=-dynamic" to the
  configure command.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              CFLAGS="" CXXFLAGS="" FFLAGS="" FCFLAGS="" \
              --enable-shared LDFLAGS=-dynamic
  ```
  See further information about "Shared and Dynamic Libraries" from the
  following NERSC web site:
  http://www.nersc.gov/users/computational-systems/edison/running-jobs/shared-and-dynamic-libraries/

* Note on compiling utility programs ncoffsets and cdfdiff.
  Error messages as shown below may appear, when the environment is set to use
  Intel-based compilers (i.e. PrgEnv-intel, default on Cray XC40).
  ```console
  /usr/bin/gcc -o ncoffsets src/utils/ncoffsets/ncoffsets.c
  In file included from /usr/include/inttypes.h:27:0,
                   from src/utils/ncoffsets/ncoffsets.c:16:
  /theta-archive/intel/compilers_and_libraries_2019.5.281/linux/compiler/include/stdint.h:43:54: error: missing binary operator before token "("
       defined(__has_include_next) && __has_include_next(<stdint.h>)
                                                        ^
  ```
  There are two ways to resolve this problem.
  1. Add "SEQ_CC=icc" to your make command line, i.e.
     make SEQ_CC=icc
  2. Unset environment variable CPATH whose value conflicts with the gcc
     include directories, i.e. run command:
     make SEQ_CC=gcc CPATH=

* Note on using compile command-line option "-fast".
  According to the NERSC URL below, the flag "-no-ipo" must be used together
  with flag "-fast" to build an application.
  http://www.nersc.gov/users/software/compilers/intel-fortran-c-and-c/

  In addition, option "-fast" implies "-static". Thus please use "-O2" flags
  when building PnetCDF shared libraries.

* Note on running "make check" or "make ptest".
  Cori is a cross-compile environment. Jobs must be run through a batch
  scheduler. At the time of this writing, the schedule is Slurm. Below
  configure command line adds 3 variables in preparation for job submission.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              CFLAGS="" CXXFLAGS="" FFLAGS="" FCFLAGS="" \
              TESTMPIRUN="srun -n NP" TESTSEQRUN="srun -n 1" \
              TESTOUTDIR="$SCRATCH"
  ```

  A sample batch script shown below allocates 4 compute nodes. "make ptest"
  requires at least 4 MPI processes. The number can be changed to use less
  number of nodes. Variable TESTOUTDIR should be set to the file directory
  for storing the test output files.
  ```console
  #!/bin/bash -l
  #SBATCH -p debug
  #SBATCH -N 4
  #SBATCH -t 00:20:00
  #SBATCH -L SCRATCH
  #SBATCH -C haswell
  make ptest
  ```

* Note on using GNU compilers
  To swap to the GNU programming environment from default, run command below
  first:
  ```console
  module swap PrgEnv-intel PrgEnv-gnu
  ```

  We found that libtool does not work well with Cray's GNU-based compiler
  wrappers, e.g. cc, CC, and ftn when building static-only library. If you
  must use GNU compilers, please add "LDFLAGS=-all-static" to the make
  command line, i.e.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              CFLAGS="" CXXFLAGS="" FFLAGS="" FCFLAGS=""
  ```

  Then add the following LDFLAGS to the make command line, i.e.
  ```console
  make       LDFLAGS=-all-static
  make tests LDFLAGS=-all-static
  make check LDFLAGS=-all-static
  ```

  Without adding LDFLAGS, you may encounter the following error messages.
  ```console
    libtool: warning: '/opt/gcc/6.3.0/snos/lib/gcc/x86_64-suse-linux/6.3.0/../../../../lib64/libgfortran.la' seems to be moved
    libtool: warning: '/opt/gcc/6.3.0/snos/lib/gcc/x86_64-suse-linux/6.3.0/../../../../lib64/libquadmath.la' seems to be moved
    libtool: warning: '/opt/gcc/6.3.0/snos/lib/gcc/x86_64-suse-linux/6.3.0/../../../../lib64/libstdc++.la' seems to be moved
    /usr/bin/ld: attempted static link of dynamic object `/opt/gcc/6.3.0/snos/lib/../lib64/libgfortran.so'
    collect2: error: ld returned 1 exit status
  ```

  * Building shared libraries:
  For 1.9.1 and later, follow the same instructions for the Intel-based
  compilers above. For 1.9.0, use command below.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc MPICXX=CC MPIF77=ftn MPIF90=ftn \
              CFLAGS="" CXXFLAGS="" FFLAGS="" FCFLAGS="" \
              --enable-shared LDFLAGS=-dynamic
  ```
  and then run
  ```console
  make
  ```
* Note on using Cray compilers (i.e. craycc, crayCC, crayftn)
  libtool fails to work with Cray compiler based wrappers. See discussion in
  https://trac.mpich.org/projects/mpich/ticket/1909 and
  http://lists.gnu.org/archive/html/libtool/2015-04/msg00003.html
  Before libtool fixes this, PnetCDF 1.9.0 does not support Cray-based
  compilers.

### Cray XC30

Building PnetCDF on the Cray XC30 (tested on Edison @ NERSC)
http://www.nersc.gov/systems/edison-cray-xc30/

The building recipe for Edison is the same as Cori.

### Cray XE6

Building PnetCDF on the Cray XE6 (tested on Hopper @ NERSC)
http://www.nersc.gov/systems/hopper-cray-xe6/

```console
./configure --prefix=/path/to/install \
            --with-mpi=/path/to/mpi/implementation \
            CFLAGS=-fast CXXFLAGS=-fast FFLAGS=-fast FCFLAGS=-fast
```

The configure command above works for PGI, GNU, and Intel
compilers, i.e. when one of the module load commands below is used:
```console
module load PrgEnv-pgi
module load PrgEnv-gnu
module load PrgEnv-intel
```
For Pathscale compilers, use commands below:
```console
module load PrgEnv-pathscale
./configure --prefix=/path/to/install \
            --with-mpi=/path/to/mpi/implementation \
            CFLAGS=-Ofast CXXFLAGS=-Ofast FFLAGS=-Ofast FCFLAGS=-Ofast
```

For Cray compilers, use commands below:
```console
module load PrgEnv-cray
./configure --prefix=/path/to/install \
            --with-mpi=/path/to/mpi/implementation \
            CFLAGS=-O2 CXXFLAGS=-O2 FFLAGS=-O2 FCFLAGS="-O2 -emf" \
            LDFLAGS=-Wl,-z,muldefs
```
Check crayftn man page for using option "-emf" in FCFLAGS: to creates .mod
files to hold module and allows the creation of lower-case module .mod file
names.

Option "-Wl,-z,muldefs" in LDFLAGS is to get around the error of multiple
definitions of `tc_version', etc.

### Cray X1

2 May 2005

I (Rob Latham) performed the following steps to get PnetCDF to build on the
Cray X1 at Oak Ridge (phoenix.ccs.ornl.gov).  Note that out-of-tree (or VPATH)
builds do not work for the Fortran interface as of 1.0.0-pre2, but we will try
to address this issue in a future release.

```console
prompt:$ module load mpt
prompt:$ export CC=cc
prompt:$ export FC=ftn
prompt:$ export MPIF77=$FC
prompt:$ export MPICC=$CC
prompt:$ export FFLAGS="-eh"
prompt:$ ./configure --prefix=/path/to/install
# note: configure takes a fairly long time.
prompt:$ make
```

The "nc_test" test will exhaust the available MPI datatypes on the X1.  Your
application might see this error:

```console
MPI has run out of internal datatype entries.
Please set the environment variable MPI_TYPE_MAX for additional space.
The current value of MPI_TYPE_MAX is 2098
```
I did as asked and nc_test completed with MPI_TYPE_MAX set to 4096

If you run on the login node, expect to see a lot of these messages:

```console
Process [nc_test] 89345 generated trap, but has signal 8 held or ignored
      epc 0x1219bb4 ra 0x1219b94 badvaddr 0x40004f0004000020
```

The messages don't *appear* to impact the program results, and additionally do
not show up if you submit the job to PBS.

Fortran codes should use '-eh' so that the Cray ftn compiler will use 1 byte
for int*1 and 2 bytes for int*2.  Otherwise, our Fortran bindings will pass
incorrect values to the C routines.

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

