#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff


# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for j in ${safe_modes} ; do
    export PNETCDF_SAFE_MODE=$j
    # echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${TESTSEQRUN} $1 arrays.bp
    ${TESTSEQRUN} $1 attributes.bp
    ${TESTSEQRUN} $1 arrays_big.bp
    echo ${ADIOS_VER_GE_1132}
    if test ${ADIOS_VER_GE_1132} = 1 ; then
       ${TESTSEQRUN} $1 attributes_big.bp
    fi
    # echo ""
done

