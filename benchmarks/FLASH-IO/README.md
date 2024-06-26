## PnetCDF FLASH-IO Benchmark

This software benchmarks the performance of PnetCDF method for the I/O and data
partitioning pattern from the FLASH I/O benchmark, developed at the Flash
Center for Computational Science, University of Chicago
(https://astro.uchicago.edu/research/flash.php).
The URL of the FLASH I/O benchmark suite is available from the following URL.
http://www.ucolick.org/~zingale/flash_benchmark_io/

The FLASH I/O benchmark suite is the I/O kernel of a block-structured adaptive
mesh hydrodynamics code developed mainly for the study of nuclear flashes on
neutron stars and white dwarfs.  The computational domain is divided into
blocks that are distributed across a number of MPI processes.  A block is a
three-dimensional array with an additional 4 elements as guard cells in each
dimension on both sides to hold information from its neighbors.  There are 24
variables per array element, and about 80 blocks on each MPI process.  A
variation in block numbers, 80, 81, or 82, per MPI process is used to simulate
a slightly unbalanced I/O load in real runs.  Since the number of blocks is
fixed for each process, increase in the number of MPI processes increases the
aggregate I/O amount as well.  FLASH I/O produces a checkpoint file containing
all 24 variables and two visualization files containing centered and corner
data covering only four plot variables.  The largest file is the checkpoint
and the I/O time of which dominates the entire benchmark.

* For some detailed description and an illustration of data partitioning pattern,
  please refer to the following paper.
  + Wei-keng Liao and Alok Choudhary. Dynamically Adapting File Domain
    Partitioning Methods for Collective I/O Based on Underlying Parallel File
    System Locking Protocols. In the Proceedings of International Conference
    for High Performance Computing, Networking, Storage and Analysis, Austin,
    Texas, November 2008.

* To compile:
  ```console
  autoreconf -i
  ./configure --with-pnetcdf=/path/PnetCDF MPIF90=/path/mpi/F90/compiler
  make
  ```

* To run:
  ```console
  mpiexec -n 512 flash_benchmark_io
  ```
  (this will use the default file base name prefix: "")
  or
  ```console
  mpiexec -n 512 flash_benchmark_io /pvfs2/flash_io_test_
  ```
  (this will use file base name prefix: "/pvfs2/flash_io_test_")

  + Command-line options:
    * [-h] print this message
    * [-q] quiet mode
    * [-b] use PnetCDF blocking APIs (default: nonblocking)
    * [-i] use MPI independent I/O (default: collective)
    * -f prefix: output file prefix name (required)


* Example output on screen:
```c
 number of guards      :             4
 number of blocks      :            80
 number of variables   :            24
 checkpoint time       :             7.52  sec
        max header     :             1.76  sec
        max unknown    :             5.76  sec
        max close      :             0.00  sec
        I/O amount     :         31109.74  MiB
 plot no corner        :             1.28  sec
        max header     :             0.29  sec
        max unknown    :             0.99  sec
        max close      :             0.00  sec
        I/O amount     :          5190.36  MiB
 plot    corner        :             1.56  sec
        max header     :             0.30  sec
        max unknown    :             1.27  sec
        max close      :             0.00  sec
        I/O amount     :          6224.35  MiB
 -------------------------------------------------------
 File base name        : /scratch2/scratchdirs/wkliao/FS_1M_128/flash_io_test_
   file striping count :           128
   file striping size  :       1048576     bytes
 Total I/O amount      :         42524.45  MiB
 -------------------------------------------------------
 nproc    array size      exec (sec)   bandwidth (MiB/s)
  512    16 x  16 x  16     10.36     4103.52
```

* Customize the problem size:
  + To change the local 3D array size, one can modify the following constants
    nxb, nyb, and nzb defined in header file physicaldata.fh.
    ```
        parameter(nxb=16)    !<<< USER EDIT
        parameter(nyb=16)    !<<< USER EDIT
        parameter(nzb=16)    !<<< USER EDIT
    ```

  + To change the maximum number of blocks each process writes, one can modify
    the variable MAXBLOCKS in file Makefile.am.
    ```
        AM_FCFLAGS += $(FC_DEFINE)MAXBLOCKS=100
    ```


*  Copyright (C) 2013, Northwestern University
*  See COPYRIGHT notice in top-level directory.

