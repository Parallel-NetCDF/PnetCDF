## Build instructions for machines available at NERSC

### [Perlmutter @ NERSC](https://docs.nersc.gov/systems/perlmutter/architecture/)

+ Below shows typical configure and make commands for building PnetCDF and its
  installation, which will produces both static and shared libraries.
  ```console
  ./configure --prefix=/path/to/install \
              cross_compiling="yes" \
              MPICC=cc \
              MPICXX=CC \
              MPIF77=ftn \
              MPIF90=ftn
  ```

+ At the time of this writing (June 2026), the default MPI compilers installed
  on Perlmutter are GNU based compilers.
  ```
  % cc --version
  gcc-14 (SUSE Linux) 14.3.0
  Copyright (C) 2024 Free Software Foundation, Inc.
  This is free software; see the source for copying conditions.  There is NO
  warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  ```

+ The Intel oneAPI compilers may not work properly when used to build the
  PnetCDF library. Note on Permutter, the GNU compiler module is loaded by
  default. If this is not your case, commands below can be used
  to switch to the GNU/NVIDIA module.
  ```console
  module swap PrgEnv-intel PrgEnv-gnu
  ```
  or
  ```console
  module swap PrgEnv-intel PrgEnv-nvidia
  ```

+ When using Intel oneAPI compilers to building PnetCDF with Fortran feature
  enabled, the building process, i.e. running "make" and "make install".
  However, using the library to compile a C application program may fail with
  the following messages complaining some MPI functions cannot be found.
  ```
  cc test.c -L /install/path/lib -lpnetcdf

  /usr/lib64/gcc/x86_64-suse-linux/14/../../../../x86_64-suse-linux/bin/ld: /opt/cray/pe/lib64/libmpifort_intel.so.12: undefined reference to `MPL_trfree'
  ```
  * In this case, user must add '-lmpifort' to the compile/link command
    line, even if the application program is not a Fortran program, for
    example,
    ```
    cc test.c -L /install/path/lib -lpnetcdf -lmpifort
    ```
  * Because of this issue, we recommend to use the GNU based compilers to
    build PnetCDF on Perlmutter.
  * This issue does not happen when the PnetCDF's Fortran feature is
    disabled, i.e. adding '--disable-fortran' at the configure command line.

---
### Note on testing PnetCDF, i.e. running "make check" or "make ptest".

Perlmutter is in a cross-compile environment, i.e. programs are compiled on the
login nodes and run on compute nodes (the hardware and software may not be the
same between the two nodes). The executable compiled for compute nodes must be
through submitting a job to a batch scheduler. At the time of this writing
(June 2026), the schedule is Slurm.

There are two `make` command options for running the test. One is "make check"
to run test programs sequentially and "make ptest" for parallel runs.

+ When running test programs developed for running on a single MPI process,
  use below lines to your batch script file which allocates a single compute
  node and appends two environment variables. Variable `TESTSEQRUN` sets the
  `srun` command-line option to use 1 process only and variable `TESTOUTDIR
  set to the file directory for storing the test output files.
  ```console
  #SBATCH -N 1
  #SBATCH -t 01:00:00

  export LD_LIBRARY_PATH="/path/to/pnetcdf/installation:$LD_LIBRARY_PATH"
  make check TESTSEQRUN="srun -n 1" \
             TESTOUTDIR="$SCRATCH"
  ```
+ When run test programs developed for running in parallel, use below lines to
  your batch script file which allocates 4 compute nodes (a less number of
  nodes is also fine).  "make ptest" requires at least 4 MPI processes. "make
  ptests" will run the tests using 2, 4, 6, and 8 processes. Variable
  `TESTMPIRUN` sets the `srun` command-line option to use `NP` MPI processes,
  where `NP` will be adjusted internally by PnetCDF. Similarly, variable
  TESTOUTDIR set to the file directory for storing the test output files.
  ```console
  #SBATCH -N 4
  #SBATCH -t 00:30:00
  make ptest TESTMPIRUN="srun -n NP" \
             TESTOUTDIR="$SCRATCH"
  ```

### Other machines historically available at NERSC
Please see [README.CRAY.md](./README.CRAY.md).

