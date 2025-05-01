#!/bin/sh
#
# Copyright (C) 2024, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# echo all commands (used for debugging)
# set -x

# This shell script requires one or two command-line options: prefx_path and
# optinal destdir_path
if [ "x$1" = x ] ; then
   echo "Usage: $0 prefx_path [destdir_path]"
   exit 1
fi

# Note command "make install prefix=/path/of/install DESTDIR=/path/to/dest"
# does not install the library. It copies all install files into folder
# "$DESTDIR/$prefix", so one can cd to $DESTDIR and pack the folder $prefix
# there into a tar ball.

prefx_path=$1
destdir_path=$2

# remove trailing '/' character
prefx_path=$(echo "$prefx_path" | sed 's:/*$::')
if [ "x$destdir_path" != x ] ; then
   destdir_path=$(echo "$destdir_path" | sed 's:/*$::')
fi

installation_path=$destdir_path$prefx_path

# echo "prefx_path=$prefx_path"
# echo "destdir_path=$destdir_path"
# echo "installation_path=$installation_path"

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
   if [ $prefixdir != $prefx_path ] ; then
       echo "Error: expecting '$prefx_path' from 'pnetcdf-config --prefix' but got $prefixdir"
       exit 1
   fi
   # check if --libdir is correctly reflecting the install path
   libdir=`$installation_path/bin/pnetcdf-config --libdir`
   if [ $libdir != $prefx_path/lib ] ; then
       echo "Error: expecting '$prefx_path/lib' from 'pnetcdf-config --libdir' but got $libdir"
       exit 1
   fi
   # check if --includedir is correctly reflecting the install path
   incdir=`$installation_path/bin/pnetcdf-config --includedir`
   if [ $incdir != $prefx_path/include ] ; then
       echo "Error: expecting '$prefx_path/include' from 'pnetcdf-config --includedir' but got $incdir"
       exit 1
   fi
fi

