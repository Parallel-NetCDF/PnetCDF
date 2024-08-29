The programs in this directory are designed to measure the I/O performance for
various of APIs as well as access patterns.

* C/aggregation.c
   + This program writes a series of 2D variables with data partitioning
     patterns of `block-block`, `*-cyclic`, `block-*`, and `*-block`,
     round-robinly.
     The `block-*` partitioning case writes 1st half followed by 2nd half. The
     same partitioning patterns are used for read. In both cases, nonblocking
     APIs are used to evaluate the performance.
   + Parameters:
     * `NVARS`: a defined C macro, the number of variables
     * `len`:   dimension size of local variables, len x len
   + Write and read performance are measured and reported separately.

* C/write_block_read_column.c
   + This program writes a series of 2D variables partitioned in a `block-block`
     pattern into a new file. The file is re-opened to read all the 2D variables
     but in a 2D `*-block` pattern.
   + Write and read performance are measured and reported separately.

* C/pnetcdf_put_vara.c
   + This program writes a series of 3D variables with 2D block-block
     partitioning pattern. Each variable is a record variable. The number of
     variables, variable size, number of time records, and the NetCDF file
     format can be customized through command-line options. In addition, option
     '-i', if set, PnetCDF nonblocking APIs will be used to write to the file.

* C/netcdf_put_vara.c
   + This sequential NetCDF-C program writes a series of 3D variables.  Each
     variable is a record variable. The number of variables, variable size,
     number of time records, and the NetCDF file format can be customized
     through command-line options.
   + This program and `C/pnetcdf_put_vara.c` can be used to compare the
     performance of NetCDF and PnetCDF when running sequentially, i.e. one
     process.

* FLASH-IO
   + FLASH is a reacting hydrodynamics code developed at University of Chicago.
     https://astro.uchicago.edu/research/flash.php
   + This benchmark is algorithmically identical to its I/O kernel.
   + This distribution contains only PnetCDF I/O method.

Copyright (C) 2012, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

