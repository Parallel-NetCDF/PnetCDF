## I/O Benchmarks
Programs in this directory are developed for performance evaluation of PnetCDF
using various I/O access patterns.

- [WRF-IO](#wrf-io)
- [FLASH-IO](#flash-io)
- [C/aggregation.c](#caggregationc)
- [C/write_block_read_column.c](#cwrite_block_read_columnc)
- [C/pnetcdf_put_vara.c](#cpnetcdf_put_varac)
- [C/netcdf_put_vara.c](#cnetcdf_put_varac)

---

### E3SM-IO
   + [E3SM](https://github.com/E3SM-Project/E3SM) (Energy Exascale Earth System
     Model) is a climate simulation model and one of the Department of Energy
     (DOE) mission applications designed to run on the DOE leadership parallel
     computers.
   + The I/O kernel of E3SM is extracted into an I/O case study in
     [E3SM-IO](https://github.com/Parallel-NetCDF/E3SM-IO) which can be used to
     evaluate the I/O performance of E3SM.
   + It's data partitioning method is based on the Hilbert space curve
     algorithm, which results in a highly irregular, non-contiguous pattern.
     Such patterns always post a great challenge for I/O libraries to achieve a
     scalable performance when storing data in files in the canonical order.
   + Example of performance evaluation results (measured on June 23, 2026).
     * Machine: [Perlmutter](https://docs.nersc.gov/systems/perlmutter/architecture/)
       at NERSC.
     * Lustre file striping configuration: 1 MB striping size and the number
       of striping count equal to the number of compute nodes allocated to
       the jobs.
     * The number of MPI processes allocated to each compute node = 128.
     * Table below shows the write bandwidths in GiB/sec in a strong scaling
       configuration. It compares PnetCDF versions 1.14.1 and 1.15.0.

       | Case   | # processes | amount (GiB) | 1.14.1 (GiB/s) | 1.15.0 (GiB/s) |
       |:------:|:-----------:|:------------:|:--------------:|:--------------:|
       | I case | 1344        |       30.0   |         0.92   |         7.16   |
       | G case | 9600        |       79.7   |         2.70   |        21.93   |
       | F case | 21600       |       28.2   |         1.19   |        12.98   |

### WRF-IO
   + [WRF](https://github.com/wrf-model/WRF) (Weather Research and Forecast
     Model) is a weather prediction computer simulation program developed at
     NCAR.
   + This benchmark is an extraction of the I/O kernel of WRF and can be used
     to evaluate the file I/O performance of WRF.
   + The metadata (dimension size, number of variables, etc.) can be found in
     [WRF-IO/wrf_header.txt)(WRF-IO/wrf_header.txt), which is also used as an
     an input file to the benchmark.
   + It's data partitioning pattern is a 2D block-block checkerboard pattern,
     along the longitude and latitude.
   + Example of performance evaluation results (measured on June 19, 2026).
     * Machine: [Perlmutter](https://docs.nersc.gov/systems/perlmutter/architecture/)
       at NERSC.
     * Lustre file striping configuration: 1 MB striping size and the number
       of striping count equal to the number of compute nodes allocated to
       the jobs.
     * The number of MPI processes allocated to each compute node = 128.
     * The number of time steps = 1
     * Table below shows the write bandwidths in GiB/sec in a weak scaling
       configuration. It compares PnetCDF versions 1.14.1 and 1.15.0.

       | # processes | grid size    | amount (GiB) | 1.14.1 (GiB/s) | 1.15.0 (GiB/s) |
       |:-----------:|:------------:|:------------:|:--------------:|:--------------:|
       | 1024        | 2600 x 3800  |       31.5   |         3.70   |         9.42   |
       | 2048        | 2600 x 7600  |       63.1   |         7.62   |        18.44   |
       | 4096        | 5200 x 7600  |      126.2   |         8.26   |        28.92   |
       | 8192        | 5200 x 15200 |      252.3   |        12.97   |        50.68   |

### FLASH-IO
   + FLASH is a reacting hydrodynamics code developed at University of Chicago.
     https://astro.uchicago.edu/research/flash.php
   + This benchmark is algorithmically identical to its I/O kernel.
   + This distribution contains only PnetCDF I/O method.

### C/aggregation.c
   + This program writes a series of 2D variables with data partitioning
     patterns of `block-block`, `*-cyclic`, `block-*`, and `*-block`,
     in a round-robin fashion.
     The `block-*` partitioning case writes 1st half followed by 2nd half. The
     same partitioning patterns are used for read. In both cases, nonblocking
     APIs are used to evaluate the performance.
   + Parameters:
     * `NVARS`: a defined C macro, the number of variables
     * `len`:   dimension size of local variables, len x len
   + Write and read performance are measured and reported separately.

### C/write_block_read_column.c
   + This program writes a series of 2D variables partitioned in a `block-block`
     pattern into a new file. The file is re-opened to read all the 2D variables
     but in a 2D `*-block` pattern.
   + Write and read performance are measured and reported separately.

### C/pnetcdf_put_vara.c
   + This program writes a series of 3D variables with 2D block-block
     partitioning pattern. Each variable is a record variable. The number of
     variables, variable size, number of time records, and the NetCDF file
     format can be customized through command-line options. In addition, option
     '-i', if set, PnetCDF nonblocking APIs will be used to write to the file.

### C/netcdf_put_vara.c
   + This sequential NetCDF-C program writes a series of 3D variables.  Each
     variable is a record variable. The number of variables, variable size,
     number of time records, and the NetCDF file format can be customized
     through command-line options.
   + This program and `C/pnetcdf_put_vara.c` can be used to compare the
     performance of NetCDF and PnetCDF when running sequentially, i.e. one
     process.

### C/indep_data_obj_create.c
   + This program create a large volume of data objects and writes a large amount of metadata
     to output file.  Each process is evenly assigned a subset of non-shared data objects
     (variables and dimensions) it intends to create.
   + The program launches MPI communications to synchronize data objects across all processes.
     All processes share a copy of the metadata so they can collectively define variables, dimensions,
     and attributes. The metadata is written to the header of output file.


Copyright (C) 2012, Northwestern University and Argonne National Laboratory.
See COPYRIGHT notice in top-level directory.

