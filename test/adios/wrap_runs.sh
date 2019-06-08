#!/bin/sh
#
# Copyright (C) 2019, Northwestern University and Argonne National Laboratory
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
    if test "$1" = ./open ; then
       ${TESTSEQRUN} $1 ${srcdir}/arrays.bp
       ${TESTSEQRUN} $1 ${srcdir}/attributes.bp
       ${TESTSEQRUN} $1 ${srcdir}/arrays_big.bp
       echo ${ADIOS_VER_GE_1132}
       if test ${ADIOS_VER_GE_1132} = 1 ; then
          ${TESTSEQRUN} $1 ${srcdir}/attributes_big.bp
       fi
    elif test "$1" = ./att ; then
       ${TESTSEQRUN} $1 ${srcdir}/attributes.bp
    else
       ${TESTSEQRUN} $1 ${srcdir}/arrays.bp
    fi
done

