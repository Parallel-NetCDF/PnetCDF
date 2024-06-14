## PnetCDF source code development repository

PnetCDF is a parallel I/O library for accessing
[Unidata's NetCDF](http://www.unidata.ucar.edu/software/netcdf) files in
classic formats. The software development is a collaborative work of
Northwestern University and Argonne National Laboratory.

* The original repository: https://svn.mcs.anl.gov/repos/parallel-netcdf
* Since June 1, 2018, PnetCDF repository has been migrated to
  here, https://github.com/Parallel-NetCDF/PnetCDF, on github.com.

### PnetCDF project web page
* https://parallel-netcdf.github.io

### PnetCDF official software releases
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
* Use the source codes in this repository only if you are interested in
  contributing the project and we welcome and appreciate any contribution.
* To clone this repository, run command
  ```console
  git clone https://github.com/Parallel-NetCDF/PnetCDF.git
  ```
  This will create a new folder named `PnetCDF`.
* Before running configure command to build PnetCDF, please run commands below
  first.
  ```console
  cd PnetCDF
  autoreconf -i
  ```
  Several files, e.g. configure and Makefile.in, will be created for running the
  configure command. The minimum versions of GNU autotools required are:
  + [autoconfig](https://www.gnu.org/software/autoconf/autoconf.html) version 2.70
  + [automake](https://www.gnu.org/software/automake) version 1.16.5
  + [libtool](https://www.gnu.org/software/libtool) version 2.4.6
  + [m4](https://www.gnu.org/software/m4/m4.html) version 1.4.17

### Build instructions and recipes
* Please read file
  [INSTALL](https://github.com/Parallel-NetCDF/PnetCDF/blob/master/INSTALL) for
  build instructions. There are also several build recipes under folder
  [doc](https://github.com/Parallel-NetCDF/PnetCDF/tree/master/doc) for a few
  popular systems.

### Current build status
* Github Actions: [![MPICH](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_mpich.yml/badge.svg)](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_mpich.yml)
[![OpenMPI](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_openmpi.yml/badge.svg)](https://github.com/Parallel-NetCDF/PnetCDF/actions/workflows/ubuntu_openmpi.yml)

### PnetCDF user guides
* [C API references](http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c)
* [Questions & Answers](http://cucis.ece.northwestern.edu/projects/PnetCDF/faq.html)
* [NetCDF4 vs. PnetCDF](./doc/netcdf4_vs_pnetcdf.md)
* [PnetCDF blocking vs. non-blocking APIs](./doc/blocking_vs_nonblocking.md)

### Mailing List
* parallel-netcdf@mcs.anl.gov
  + To subscribe, please visit the list information page
    https://lists.mcs.anl.gov/mailman/listinfo/parallel-netcdf
  + The past discussions in the mailing list are available in:
    http://lists.mcs.anl.gov/pipermail/parallel-netcdf/.

### Project funding supports:
* Ongoing development and maintenance of PnetCDF is supported by the [Exascale
  Computing Project](https://www.exascaleproject.org) (17-SC-20-SC), a joint
  project of the U.S. Department of Energy's Office of Science and National
  Nuclear Security Administration, responsible for delivering a capable
  exascale ecosystem, including software, applications, and hardware
  technology, to support the nation's exascale computing imperative.

