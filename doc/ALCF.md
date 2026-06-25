## Build instructions for machines available at ALCF

### [Aurora @ ALCF](https://docs.alcf.anl.gov/aurora)

* Below shows typical configure and make commands for building PnetCDF and its
  installation, which will produces both static and shared libraries.
  ```console
  ./configure --prefix=/path/to/install \
              MPICC=mpicc \
              MPICXX=mpicxx \
              MPIF77=mpifort \
              MPIF90=mpifort \
  ```
* At the time of this writing (June 2026), the default MPI compilers installed
  on Aurora are Intel oneAPI based compilers.
  ```
  % mpicc --version
  Intel(R) oneAPI DPC++/C++ Compiler 2025.3.2 (2025.3.2.20260112)
  Target: x86_64-unknown-linux-gnu
  Thread model: posix
  InstalledDir: /opt/aurora/26.26.0/oneapi/compiler/2025.3/bin/compiler
  Configuration file: /opt/aurora/26.26.0/oneapi/compiler/2025.3/bin/compiler/../icx.cfg
  ```

---
### Note on testing PnetCDF, i.e. running "make check" or "make ptest".

Aurora is in a cross-compile environment, i.e. programs are compiled on the
login nodes and run on compute nodes (the hardware and software may not be the
same between the two nodes). The executable compiled for compute nodes must be
through submitting a job to a batch scheduler. At the time of this writing
(June 2026), the batch job scheduler is PBS (Portable Batch System).

There are two `make` command options for running the test. One is "make check"
to run test programs sequentially and "make ptest" for parallel runs.

+ When running test programs developed for running on a single MPI process,
  use below lines to your batch script file which allocates a single compute
  node and appends two environment variables. Variable `TESTSEQRUN` sets the
  `mpiexec` command-line option to use 1 process only and variable `TESTOUTDIR
  set to the file directory for storing the test output files.
  ```console
  #PBS -l select=1
  #PBS -l walltime=01:00:00

  export LD_LIBRARY_PATH="/path/to/pnetcdf/installation:$LD_LIBRARY_PATH"

  make check TESTSEQRUN="mpiexec -n 1" \
             TESTOUTDIR="/output/folder"
  ```
+ When run test programs developed for running in parallel, use below lines to
  your batch script file which allocates 4 compute nodes (a less number of
  nodes is also fine).  "make ptest" requires at least 4 MPI processes. "make
  ptests" will run the tests using 2, 4, 6, and 8 processes. Variable
  `TESTMPIRUN` sets the `mpiexec` command-line option to use `NP` MPI
  processes, where `NP` will be adjusted internally by PnetCDF. Similarly,
  variable TESTOUTDIR set to the file directory for storing the test output
  files.
  ```console
  #SBATCH -N 4
  #SBATCH -t 00:30:00
  make ptest TESTMPIRUN="mpiexec -n NP" \
             TESTOUTDIR="/output/folder"
  ```

---
### [Polaris @ ALCF](https://docs.alcf.anl.gov/polaris)

* Below shows typical configure and make commands for building PnetCDF and its
  installation, which will produces both static and shared libraries.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc \
              MPICXX=CC \
              MPIF77=ftn \
              MPIF90=ftn
  ```

* At the time of this writing (June 2026), the default MPI compilers installed
  on Polaris are NVIDIA based compilers.
  ```
  % cc --version

  nvc 25.5-0 64-bit target on x86-64 Linux -tp zen3-64
  NVIDIA Compilers and Tools
  Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
  ```

* The batch job scheduler on Polaris is also PBS (Portable Batch System), same
  as on Aurora. Please refer to the above note for running test on Polaris.

### Other machines historically available at ALCF
Please see [README.CRAY.md](./README.CRAY.md).

