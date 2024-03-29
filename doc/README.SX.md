## NEC SX

* Current notes for NEC SX; based on pnetcdf version 1.0.0; 28 July, 2005

* SX Cross compiler environment is not supported yet (the check for a working
  ftruncate is not possible). Configure steps have to be invoked on the SX directly.
  Configuring on the target host and then cross-compiling works fine when setting
  up necessary aliases on the configure host (cc -> sxcc, etc) or vice versa
  on the compile host (sxcc -> cc).

* With the following environment variables a libpnetcdf.a has been built
  successfully
  ```console
    MPICC=mpic++
    MPIF77=mpif90
    FC=f90
    CC=c++
    FFLAGS=-dW
  ```

* Built on NEC SX6 with:
   - Operating system SUPER-UX 14.1
   - C++/SX compiler rev.061 2004/01/06
   - f90/SX compiler rev.285 2003/09/25

* Note on nf_test
  + `-dW` disables promotion of Integer*2. However no interfaces for
    `Integer*1` are built. This causes 2 compile time errors in util.F. Lines
    1158 and 1226 have to be turned into valid Fortran syntax.

* Note on test_dtype
  + The executables have to be compiled with `-pvctl loopcnt=186048`. To avoid
    run time error in test_array the vectorisation of the loop starting at line
    256 needs to be disabled with `#pragma cdir novector`

* Rob Latham <robl@mcs.anl.gov>
* Rene Redler <redler@ccrl-nece.de>
* Joachim Worringen <joachim@ccrl-nece.de>

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.
