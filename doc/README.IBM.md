## Build PnetCDF using the IBM XL Compilers

1. Building PnetCDF on BGQ
2. Building PnetCDF on BGP
3. Building PnetCDF on BGL
4. Building PnetCDF on UCAR BlueFire
5. Building PnetCDF on IBM SP

### Building PnetCDF on BGQ

Building for BGQ is not so different from BGP or BGL: front end node is still a
cross compile host for the back end (as of this writing, the "log in io nodes",
a.k.a lion nodes, are not on line for the Mira BGQ machine).

Be sure to run configure with the --build and --host flags to put it in "cross
compile mode".  This will make configure use compile-only tests, instead of the
usual compile-and-run tests (running tests on the BGP login node won't work as
the compute nodes are totally different. This will change if/when lion nodes
are available).

* To build PnetCDF with the IBM xl compilers
  ```console
  ./configure --prefix=/path/to/install \
              --host powerpc64-bgq-linux \
              --build ppc64-redhat-linux \
              MPICC=mpixlc \
              MPICXX=mpixlcxx \
              MPIF77=mpixlf77 \
              MPIF90=mpixlf90 \
              FFLAGS="-qsuppress=cmpmsg" \
              FCFLAGS="-qsuppress=cmpmsg" \
  ```
Fortran flag "-qsuppress=cmpmsg" can suppresses XLF compile informational
messages that report compilation progress and a successful completion.

If you encounter the following warning message,
"(W) Option DLINES is ignored within Fortran 90 free form and IBM free form"
it is safe to ignore. This is due to Automake adds -DHAVE_CONFIG_H to Fortran
compile command line, in addition to C. However, xlf compiler uses -WF,-D
instead of -D for command-line macro definition.

If you encounter error messages related to Darshan when running "make check",
"make ptest", or "make tests", such as
```console
/soft/perftools/darshan/darshan-2.3.0/lib/libdarshan-mpi-io.a(darshan-pnetcdf.o): In function `__wrap_ncmpi_close':
try disable Darshan by setting the environment variable DARSHAN_DISABLE to 1.
```

* To build shared library with XL based compilers
  + PnetCDF's default is to build static libraries only. To build shared
    libraries, add the following to the above configure command line.
    ```console
                --enable-shared \
                LDFLAGS=-qmkshrobj
    ```
  + Link flag -qmkshrobj is required to build shared libraries when using the
    IBM xl C/C++ compilers. It should not be used when building static
    libraried only.

  + See the Mira compiling and linking FAQ on shared libraries at
    https://www.alcf.anl.gov/user-guides/compiling-and-linking-faq#how-do-i-compile-and-link-a-shared-library

* To build PnetCDF with the GNU compilers

  + At first, you may want to "swap" the softenv to use the GNU compiler
    wrappers.  The commands to do so are:
    ```console
    soft delete +mpiwrapper-xl
    soft add +mpiwrapper-gcc
    ```
  + The GNU-based MPI compilers are mpicc, mpicxx, mpif77, and mpif90
  + To build static-only library with GNU based compilers
    ```console
    ./configure --prefix=/path/to/install \
                --host powerpc64-bgq-linux \
                --build ppc64-redhat-linux \
                MPICC=mpicc \
                MPICXX=mpicxx \
                MPIF77=mpif77 \
                MPIF90=mpif90
    ```
  + Then add the following LDFLAGS to the make command line, i.e.
    ```console
        make       LDFLAGS=-all-static
        make tests LDFLAGS=-all-static
        make check LDFLAGS=-all-static
    ```
* To build shared library with GNU based compilers
  + PnetCDF's default is to build static libraries only. To build shared
    libraries, add the following to the above configure command line.
    ```console
                --enable-shared \
                LDFLAGS=-dynamic
    ```
  + Link flag -dynamic is required to build shared libraries when using the
    GNU C/C++ compilers. It should not be used when building static libraries
    only.

* To use thread-safe XL compilers
  + You can simply replace the XL based MPI compilers with the followings
    in the configure command line.
    ```console
    MPICC=mpixlc_r    \
    MPICXX=mpixlcxx_r \
    MPIF77=mpixlf77_r \
    MPIF90=mpixlf90_r \
    ```

* As far as we know, there are no issues with PnetCDF on BGQ, but should you
  find any, please email parallel-netcdf@mcs.anl.gov

* Note on compiling Fortran 90 programs on Intrepid @ ANL
https://www.alcf.anl.gov/user-guides/bgp-faqs-compiling-and-linking#how-do-i-use-xlf-and-"use-mpi"

* If you attempt to do a "use mpi" from XL Fortran, and build with the default
  compiler wrappers, you will encounter this error:
  ```console
  XL Fortran90 error: Module symbol file for module mpi is in a format not
recognized by this compiler
  ```

* The cause of the error is that the mpi.mod referenced in the "default"
  software stack (/bgsys/drivers/ppcfloor/comm/default) is built with GNU
  Fortran.
* A work-around to the problem using mpif.h instead of "use mpi.mod", e.g.:
  ```
  program test
  include "mpif.h"

  ... your MPI stuff ...

  end program test
  ```

* The proper solution is to use an alternate software stack with MPI built
  using the XL compilers instead of the GNU compilers, thus containing a
  Fortran mpi module (mpi.mod) generated by XLF. This software stack was
  introduced in V1R4M1. You may build with it using the wrappers in comm/xl:
  ```
 /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf90
  ```

### Note on running PnetCDF testing programs on Mira/Cetus @ ALCF, ANL

Mira and Cetus are IBM BlueGene Q, a cross-compile system. Hence, one cannot
run PnetCDF test programs on the login nodes. Test programs must run on the
compute nodes through submitting a batch job. An example batch script file and
the submission command for allocating 8 node with 90 minutes are shown below.
Note the testing requires at least 8 MPI processes and in this example the
command requests 8 nodes and the script runs 1 MPI process per node.

```console
% cat cobalt.script
#!/bin/sh
echo "Starting Cobalt job script"

# test seqential programs
make check TESTMPIRUN="runjob --block $COBALT_PARTNAME --ranks-per-node 1 --np NP : " \
           TESTSEQRUN="runjob --block $COBALT_PARTNAME --ranks-per-node 1 --np 1 : " \
           TESTOUTDIR=/path/to/GPFS/directory

# test parallel programs
make ptest TESTMPIRUN="runjob --block $COBALT_PARTNAME --ranks-per-node 1 --np NP : " \
           TESTSEQRUN="runjob --block $COBALT_PARTNAME --ranks-per-node 1 --np 1 : " \
           TESTOUTDIR=/path/to/GPFS/directory


% qsub -A YOUR_ACCOUNT -n 8 -t 30 --mode script ./cobalt.script
```

* Running "make testing" can take a while (more than one hour on Mira/Cetus).
  This is because the testing performs many small read/write requests and under
  BlueGene's I/O forwarding architecture, local file-system caching is not
  available on the compute nodes.

* Running "make ptest" will take approximately 10 minutes.

### Building PnetCDF on BGP

* Building for BGP is not so different from BGL: front end node is still a
  cross compile host for the back end.

* Be sure to run configure with the --build and --host flags to put it in
  "cross compile mode".  This will make configure use compile-only tests,
  instead of the usual compile-and-run tests (running tests on the bgp login
  node won't work as the compute nodes are totally different).

* There is one run-time check for MPI-IO support of resized types.
  Unfortunately we have to test for this with a runtime test, but you can set
  the environment variable "ac_cv_MPI_TYPE_RESIZED_WORKS" .
  ```console
  ./configure --host powerpc-bgp-linux --build powerpc64-suse-linux \
              --with-mpi=/bgsys/drivers/ppcfloor/comm
  ```

* It's possible to build PnetCDF with the IBM xl compilers:
  ```console
  ./configure --host powerpc-bgp-linux --build powerpc64-suse-linux \
              --with-mpi=/bgsys/drivers/ppcfloor/comm/xl \
              MPICC=mpixlc MPICXX=mpixlcxx MPIF77=mpixlf77 MPIF90=mpixlf90 \
              FFLAGS="-qsuppress=cmpmsg" FCFLAGS="-qsuppress=cmpmsg"
  ```

* Fortran flag "-qsuppress=cmpmsg" can suppresses XLF compile informational
  messages that report compilation progress and a successful completion.

* As far as we know, there are no issues with PnetCDF on BGP, but should you
  find any, email parallel-netcdf@mcs.anl.gov

### Note on compiling Fortran 90 programs on Intrepid @ ANL
https://www.alcf.anl.gov/user-guides/bgp-faqs-compiling-and-linking#how-do-i-use-xlf-and-"use-mpi"

* If you attempt to do a "use mpi" from XL Fortran, and build with the default
  compiler wrappers, you will encounter this error:
  ```console
  XL Fortran90 error: Module symbol file for module mpi is in a format not
recognized by this compiler
  ```

* The cause of the error is that the mpi.mod referenced in the "default"
  software stack (/bgsys/drivers/ppcfloor/comm/default) is built with GNU
  Fortran.
* A work-around to the problem using mpif.h instead of "use mpi.mod", e.g.:
  ```console
  program test
  include "mpif.h"

  ... your MPI stuff ...

  end program test
  ```

* The proper solution is to use an alternate software stack with MPI built
  using the XL compilers instead of the GNU compilers, thus containing a
  Fortran mpi module (mpi.mod) generated by XLF. This software stack was
  introduced in V1R4M1. You may build with it using the wrappers in comm/xl:
  ```console
 /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf90
  ```

### Building PnetCDF on BGL

Be sure to run configure with the --build and --host flags to put it in "cross
compile mode".  This will make configure use compile-only tests, instead of the
usual compile-and-run tests (running tests on the bgl login node won't work as
the compute nodes are totally different).

```console
configure --build powerpc32-unknown-gnu --host powerpc-suse-linux  \
          --with-mpi=/bgl/BlueLight/ppcfloor/bglsys/
```

* It's possible to build PnetCDF with the IBM xl compilers, but you have to set
  quite a few environment variables

  ```console
  export CC=blrts_xlc
  export MPICC=blrts_xlc
  export CXX=blrts_xlC
  export MPICXX=blrts_xlC
  export FC=blrts_xlf
  export F77=blrts_xlf
  export MPIF77=blrts_xlf

  export CFLAGS="-I/bgl/BlueLight/ppcfloor/bglsys/include"
  export LIBS="-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts -lrts.rts"

  configure --build powerpc32-unknown-gnu --host powerpc-suse-linux
  ```

* Several early versions of IBM's MPI-IO implementation would segfault under
  certain workloads.  If you are running driver version
  "V1R3M0_240_2006-060623" or newer, the segfault issue should be resolved.  If
  you are running an older driver, read on:

* When built against some older BlueGene drivers,  nc_test does not run
  completely without setting a special environment variable, hitting a seg
  fault deep inside ROMIO.  We first noticed this in IBM's "Driver 202" MPI and
  also in "V1R1M1_253_2005-051003" and "V1R2M1_020_2006-060110"  We have told
  IBM developers about the problem.  Code that makes use of the
  ncmpi_get_vara_*_all or ncmpi_put_vara_*_all routines will likely trigger a
  seg fault.  IBM has provided a workaround:  if your code seg-faults, try
  setting the "BGLMPIO_TUNEBLOCKING" environment variable to 0.  With this
  environment variable set, nc_tests runs to completion and passes.  For one
  real-world example, the FLASH-IO benchmark with 8 processors sees a 5%
  performance hit when writing out plotfiles. That increases to 22% with 16
  processors. Again, upgrading to the latest BlueGene drivers should fix this
  issue.

### Building PnetCDF on UCAR BlueFire
* the nc_test and nf_test tests use M4 to generate code.  AIX-m4 won't
  generate correct code, so use gnu M4.
* We had a report of a build failure when building the 'ncmpigen' utility.  The
  system 'bison' tool (bison-1.875) did not produce output files, but when the
  user installed a new version of gnu bison, everything worked fine.
* the PnetCDF code is slowly taking on more and more c99 features.
* On BlueFire, the commands I used look like this:
  ```console
  $ module add m4-1.4.14

  $ configure --prefix=/path/to/install \
          CFLAGS=-qlanglvl=stdc99 \
          CC=xlc FC=xlf F77=xlf F90=xlf90 \
          MPICC=mpcc_r MPIF77=mpxlf_r MPIF90=mpxlf90_r CXX=xlC MPICXX=mpCC_r
  ```

### Building PnetCDF on IBM SP
* John Tannahill <tannahill1@llnl.gov> reported success building PnetCDF on the
  'seaborg' cluster (an IBM-SP at NERSC) and Tyce Mclarty <mclarty3@llnl.gov>
  reported success on LLNL's 'frost' cluster by setting these environment
  variables:
  ```console
   setenv MPICC  mpcc_r
   setenv MPIF77 mpxlf_r
   setenv F77    xlf
   setenv FC     xlf
   setenv CC     xlc
   setenv CXX    xlC
  ```
* after setting these variables, configure/make/make install should "just
  work".  We also successfully tested 64-bit mode on 'DataStar' (an IBM-SP at
  SDSC) with environment variables:
  ```console
   setenv OBJECT_MODE 64
   setenv MPICC  mpcc_r
   setenv MPIF77 mpxlf_r
   setenv F77    xlf
   setenv FC     xlf
   setenv CC     xlc
   setenv CXX    xlC
   setenv CFLAGS -q64
   setenv FFLAGS -q64
   setenv F90FLAGS -q64
   setenv CXXFLAGS -q64
  ```

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

