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
    ${TESTSEQRUN} $1              ${TESTOUTDIR}/pres_temp_4D.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/pres_temp_4D.nc
done

if [ -n "${TESTBB}" ]; then
   for j in 0 1 ; do
       export PNETCDF_SAFE_MODE=$j
       export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/pres_temp_4D.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/pres_temp_4D.nc
   done
fi
