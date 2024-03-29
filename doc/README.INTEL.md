## Build PnetCDF Using The Intel Compilers

* When using old versions of Intel MPI compilers (4.x), you might encounter the
  following error message.
  ```console
  In file included from /opt/intel/impi/4.1.3/include/mpi.h(1279),
                   from ../lib/pnetcdf.h(10),
                   from ../../../trunk/src/libcxx/ncmpiType.h(2),
                   from ../../../trunk/src/libcxx/ncmpiType.cpp(2):
  /opt/intel/impi/4.1.3/include/mpicxx.h(95): error: #error directive: "SEEK_SET is
  #defined but must not be for the C++ binding of MPI. Include mpi.h before stdio.h"
  #error "SEEK_SET is #defined but must not be for the C++ binding of MPI. Include mpi.h before stdio.h"
     ^
  ```
* A solution is to add the following C++ preprocessor flags at the configure
  command line, for example.
  ```console
  ./configure CXXCPPFLAGS="-DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX"
  ```

* Note that this error has been resolved when using newer Intel MPI compiler
  (5.x).  See the two URLs below for more information.
  + https://software.intel.com/en-us/articles/intel-cluster-toolkit-for-linux-error-when-compiling-c-aps-using-intel-mpi-library-compilation-driver-mpiicpc
  + https://wiki.mpich.org/mpich/index.php/Frequently_Asked_Questions#Q:_I_get_compile_errors_saying_.22SEEK_SET_is_.23defined_but_must_not_be_for_the_C.2B.2B_binding_of_MPI.22

* If you encountered the following error when using mpiicc compiler, add
  "CPP="mpiicc -E"" to your configure command line should resolve the problem.
  ```console
  configure: error: in pnetcdf-1.11.0': configure: error: C preprocessor "mpiicc" fails sanity check See config.log' for more details
  ```
  + See issue [#43](https://github.com/Parallel-NetCDF/PnetCDF/issues/43).

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

