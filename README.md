## PnetCDF source code development repository
[![MPICH](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_mpich.yml/badge.svg)](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_mpich.yml)
[![OpenMPI](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_openmpi.yml/badge.svg)](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_openmpi.yml)


PnetCDF is a parallel I/O library for accessing
[Unidata's NetCDF](http://www.unidata.ucar.edu/software/netcdf) files in
classic formats. The software development is a collaborative work of
Northwestern University and Argonne National Laboratory.

* The original repository: https://svn.mcs.anl.gov/repos/parallel-netcdf
* Since June 1, 2018, the PnetCDF repository has been migrated to
  here on github.

### PnetCDF project home page
* https://parallel-netcdf.github.io
  contains more information about PnetCDF.

### PnetCDF official software releases
* The latest stable release is
  [pnetcdf-1.14.1.tar.gz](https://parallel-netcdf.github.io/Release/pnetcdf-1.14.1.tar.gz)
  ([release note](https://github.com/Parallel-NetCDF/Parallel-NetCDF.github.io/blob/master/Release_notes/1.14.1.md)),
  available since July 31, 2025.
* All **official released versions** can be found in
  https://parallel-netcdf.github.io/wiki/Download.html
* Note the ["releases"](https://github.com/Parallel-NetCDF/PnetCDF/releases)
  link on this page above contains only tagged versions. They are by no means
  official releases, but simply checkpoints. They contain unused historical
  files. Users are recommended to download the official releases, not tagged
  versions.

### Build PnetCDF using source codes of this development repository
* The source codes in this repository are constantly under development. They
  should **NOT** be used for production runs.
* Building PnetCDF using the source codes in this repository only if you are
  interested in contributing the project and we welcome and appreciate any
  contribution.
* Required software
  + MPI compilers, e.g. [MPICH](https://www.mpich.org) and
    [OpenMPI](https://www.open-mpi.org)
  + [autoconfig](https://www.gnu.org/software/autoconf) version 2.71
  + [automake](https://www.gnu.org/software/automake) version 1.16.5
  + [libtool](https://www.gnu.org/software/libtool) version 2.5.4
  + [m4](https://www.gnu.org/software/m4) version 1.4.17
* Clone, build, and installation commands:
  ```console
  git clone https://github.com/Parallel-NetCDF/PnetCDF.git

  cd PnetCDF
  autoreconf -i

  ./configure --prefix=/install/path --with-mpi=/mpi/path
  make install
  ```

### Build instructions and recipes
* Please read file [INSTALL](./INSTALL) for build instructions. There are also
  several build recipes under folder [doc](./doc#readme) for a few popular
  systems.

### Users documents
* [C API references](https://parallel-netcdf.github.io/doc/c-reference/pnetcdf-c/index.html)
* [Questions & Answers](https://parallel-netcdf.github.io/doc/faq.html)
* [Utility Programs](./src/utils#readme)
* [NetCDF4 vs. PnetCDF](./doc/netcdf4_vs_pnetcdf.md)
* PnetCDF [blocking vs. non-blocking APIs](./doc/blocking_vs_nonblocking.md)
* [CDL header API references](./doc/cdl_api_guide.md)

### Mailing List
* parallel-netcdf@mcs.anl.gov
  + To subscribe, please visit the list information page
    https://lists.mcs.anl.gov/mailman/listinfo/parallel-netcdf
  + The past discussions in the mailing list are available in:
    http://lists.mcs.anl.gov/pipermail/parallel-netcdf/.

### External links to some related projects and application users
* [PnetCDF Python](https://github.com/Parallel-NetCDF/PnetCDF-Python),
  a python interface to PnetCDF library.
* [E3SM I/O kernel study](https://github.com/Parallel-NetCDF/E3SM-IO)
* [Scorpio](https://github.com/E3SM-Project/scorpio), the I/O module of E3SM.
* [PIO](https://github.com/NCAR/ParallelIO) - parallel I/O library at NCAR.
* [WRF](https://github.com/wrf-model/WRF/tree/master/external/io_pnetcdf) -
  Weather Research & Forecasting Model at NCAR.


### Acknowledgements
* Ongoing development and maintenance of PnetCDF-python is supported by the
  U.S. Department of Energy's Office of Science, Scientific Discovery through
  Advanced Computing (SciDAC) program, OASIS Institute.
* From 2016 to 2023, the development and maintenance of PnetCDF was supported
  by the [Exascale Computing Project](https://www.exascaleproject.org)
  (17-SC-20-SC), a joint project of the U.S. Department of Energy's Office of
  Science and National Nuclear Security Administration, responsible for
  delivering a capable exascale ecosystem, including software, applications,
  and hardware technology, to support the nation's exascale computing
  imperative.
* The PnetCDF project has been continuously supported by the U.S. Department of
  Energy's Office of Science, Scientific Discovery through Advanced Computing
  (SciDAC) program since its initiation in 2001.


