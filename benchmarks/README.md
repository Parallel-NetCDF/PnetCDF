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

* FLASH-IO
   + FLASH is a reacting hydrodynamics code developed at University of Chicago.
     https://astro.uchicago.edu/research/flash.php
   + This benchmark is algorithmically identical to its I/O kernel.
   + This distribution contains only PnetCDF I/O method.

