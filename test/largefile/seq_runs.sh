#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

for j in 0 ; do
    # disable safe mode, as large tests already run slow
    export PNETCDF_SAFE_MODE=$j
    echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

    for i in ${TESTPROGRAMS}; do
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/$i.nc
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
    done
done
