## PnetCDF source code repository

PnetCDF is a parallel I/O library for accessing [NetCDF](http://www.unidata.ucar.edu/software/netcdf) files. The software development is a collaborative work of Northwestern University and Argonne National Laboratory.

### PnetCDF project web page
* https://parallel-netcdf.github.io

### PnetCDF official releases of source codes
* All official released versions can be found in http://cucis.ece.northwestern.edu/projects/PnetCDF/download.html
* Note the github link marked "releases" in this page contains only tagged versions. They are by no means official releases, but simply checkpoints. Users are recommended to use the official releases, not tagged versions.

### Build PnetCDF using source codes of this repository
* The source codes in this repository are under development. They should not be used for production runs.
* To clone this repository, run command
  ```
  git clone https://github.com/Parallel-NetCDF/PnetCDF.git
  ```
  This will create a new folder named `PnetCDF`.
* Before running configure command to build PnetCDF, please run command below first.
  ```
  cd PnetCDF
  autoreconf -i
  ```
  Several autotools files will be created to be used by configure command.
  The minimum versions of GNU autotools required are:
  + [autoconfig](https://www.gnu.org/software/autoconf/autoconf.html) version 2.69
  + [automake](https://www.gnu.org/software/automake) version 1.13
  + [libtool](https://www.gnu.org/software/libtool) version 2.4.2
  + [m4](https://www.gnu.org/software/m4/m4.html) version 1.4.17

### Build instructions and recipes
* Please read file INSTALL for build instructions. There are also several build recipes under folder `doc` for a few popular systems.

### Current build status
* [Travis CI ![Build Status](https://travis-ci.org/Parallel-NetCDF/PnetCDF.svg?branch=master)](https://travis-ci.org/Parallel-NetCDF/PnetCDF)
* [Coverity Scan ![Build Status](https://scan.coverity.com/projects/15801/badge.svg)](https://scan.coverity.com/projects/parallel-netcdf-pnetcdf)

### PnetCDF user guide
* C references: http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c
* Questions & Answers: http://cucis.ece.northwestern.edu/projects/PnetCDF/faq.html

### Mailing List
* parallel-netcdf@mcs.anl.gov
  + To subscribe, please visit the list information page https://lists.mcs.anl.gov/mailman/listinfo/parallel-netcdf
  + The past discussions in the mailing list are available in: http://lists.mcs.anl.gov/pipermail/parallel-netcdf/.

### Project funding supports:
* Ongoing development and maintenance of PnetCDF is supported by the [Exascale
  Computing Project](https://www.exascaleproject.org) (17-SC-20-SC), a joint
  project of the U.S. Department of Energy's Office of Science and National
  Nuclear Security Administration, responsible for delivering a capable
  exascale ecosystem, including software, applications, and hardware
  technology, to support the nation's exascale computing imperative.

