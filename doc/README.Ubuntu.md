## Build PnetCDF on Ubuntu

The following waning message may occur during "make". It is safe to ignore.
* See http://gnu-automake.7480.n7.nabble.com/bug-20082-new-warning-from-ar-on-rawhide-systems-td21769.html
```console
ar: `u' modifier ignored since `D' is the default (see `U')
```

### The following note applies to 1.9.0 only, when building shared library.

* Set LDFLAGS to the following at the configure command line, i.e.,
  ```console
  ./configure --prefix=/path/to/install \
              --enable-shared \
              LDFLAGS="-Wl,--allow-shlib-undefined"
  ```

* Without option "--enable-shared", only the static libraries will be built
  and in this case there is no need to set LDFLAGS.

* Note that without setting the LDFLAGS to the above, you might see error
  messages similar to below, when using gfortran based MPI compilers.
  ```console
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_shape_4'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_os_error'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_runtime_error_at'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_runtime_error'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_internal_unpack'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_compare_string'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_internal_pack'
   ../../../src/libs/.libs/libpnetcdf.so: undefined reference to `_gfortran_string_len_trim'
   collect2: error: ld returned 1 exit status
  ```

Copyright (C) 2017, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

