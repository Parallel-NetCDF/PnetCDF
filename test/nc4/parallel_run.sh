#!/bin/bash
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

PNETCDF_HINTS=
if test "x$GIO_ONLY" = x1 ; then
   PNETCDF_HINTS+="nc_driver=gio;"
fi

if test "x$MIMIC_LUSTRE" != x1 ; then
   PNETCDF_HINTS+="cb_nodes=2"
fi
export PNETCDF_HINTS=$PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc.cdf4
done

