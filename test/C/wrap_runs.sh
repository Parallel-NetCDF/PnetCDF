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
    if [ -n "${TESTNC4}" ]; then
        ${TESTSEQRUN} $1              ${TESTOUTDIR}/pres_temp_4D_nc4.nc 4
        # Validator does not support nc4
        # ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/pres_temp_4D_nc4.nc 4
    fi
done
