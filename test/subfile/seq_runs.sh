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
    for i in $TESTPROGRAMS; do
        ${TESTSEQRUN} ./$i         -f ${TESTOUTDIR}/$i.nc -s 2
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
    done
done

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test_subfile.nc.subfile_0.nc


if [ -n "${TESTBB}" ]; then
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        for i in $TESTPROGRAMS; do
            export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable"
            ${TESTSEQRUN} ./$i         -f ${TESTOUTDIR}/$i.nc -s 2
            unset PNETCDF_HINTS
            ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        done
    done

    export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable"
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test_subfile.nc.subfile_0.nc
    unset PNETCDF_HINTS
fi
