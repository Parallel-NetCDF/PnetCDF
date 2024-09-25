#!/bin/sh
#
# Copyright (C) 2024, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# This requires one command-line option, the installation path
if [ "x$1" = x ] ; then
   echo "Usage: $0 /pnetcdf/install/path"
   exit 1
fi

installation_path=$1

# check if folder installation_path exists
if [ ! -d $installation_path ]; then
    echo "Error: folder $installation_path cannot be found"
    exit 1
fi

# remove trailing '/' character
# Note on Mac OSX, realpath does not support option -s
# Ideally, -s should be used to avoid expanding a symlink
installation_path=$(realpath $installation_path)

# check if pnetcdf_version exists in the install folder
if [ ! -x $installation_path/bin/pnetcdf_version ]; then
    echo "Error: pnetcdf_version not found in $installation_path/bin"
    exit 1
fi

# check if pnetcdf.pc exists in the install folder
if [ ! -f $installation_path/lib/pkgconfig/pnetcdf.pc ]; then
    echo "Error: pnetcdf.pc not found in $installation_path/lib/pkgconfig"
    exit 1
fi

# check if pnetcdf-config exists in the install folder
if [ ! -x $installation_path/bin/pnetcdf-config ]; then
    echo "Error: pnetcdf-config not found in $installation_path/bin"
    exit 1
else
   # check if --prefix is correctly reflecting the install path
   prefixdir=`$installation_path/bin/pnetcdf-config --prefix`
   if [ $prefixdir != $installation_path ] ; then
       echo "Error: expecting '$installation_path' from 'pnetcdf-config --prefix' but got $prefixdir"
       exit 1
   fi
   # check if --libdir is correctly reflecting the install path
   libdir=`$installation_path/bin/pnetcdf-config --libdir`
   if [ $libdir != $installation_path/lib ] ; then
       echo "Error: expecting '$installation_path/lib' from 'pnetcdf-config --libdir' but got $libdir"
       exit 1
   fi
   # check if --includedir is correctly reflecting the install path
   incdir=`$installation_path/bin/pnetcdf-config --includedir`
   if [ $incdir != $installation_path/include ] ; then
       echo "Error: expecting '$installation_path/include' from 'pnetcdf-config --includedir' but got $incdir"
       exit 1
   fi
fi

