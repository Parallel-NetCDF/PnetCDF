#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

for j in 0 1 ; do
    export PNETCDF_SAFE_MODE=$j
    # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${TESTSEQRUN} $1              ${TESTOUTDIR}/$1.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$1.nc
done

if [ -n "${TESTDW}" ]; then
   for j in 0 1 ; do
       export PNETCDF_SAFE_MODE=$j
       export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/$1.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$1.nc
   done
fi

