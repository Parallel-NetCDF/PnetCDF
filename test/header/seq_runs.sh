#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# header consistency tests are designed to run on more than one MPI process
for j in 0 1 ; do
    export PNETCDF_SAFE_MODE=$j
    echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    for i in ${TESTPROGRAMS} ; do
        ${TESTSEQRUN} ./$i ${TESTOUTDIR}/$i.nc
    done
done

echo ""

if [ -n "${TESTBB}" ]; then
    echo "---- testing burst buffering"
    export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        for i in ${TESTPROGRAMS} ; do
            ${TESTSEQRUN} ./$i ${TESTOUTDIR}/$i.nc
        done
    done
    unset PNETCDF_HINTS
fi
