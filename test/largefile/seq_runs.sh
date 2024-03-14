#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

# disable safe mode, as large tests already run slow
export PNETCDF_SAFE_MODE=0

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${TESTPROGRAMS}; do
    ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/$i.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc
done
