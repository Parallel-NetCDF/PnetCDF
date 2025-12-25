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
if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for j in ${safe_modes} ; do
    if test "$j" = 1 ; then # test only in safe mode
       SAFE_HINTS="romio_no_indep_rw=true"
    else
       SAFE_HINTS="romio_no_indep_rw=false"
    fi
    PNETCDF_HINTS=
    if test "x$SAFE_HINTS" != x ; then
       PNETCDF_HINTS="$SAFE_HINTS"
    fi
    export PNETCDF_HINTS="$PNETCDF_HINTS"
    export PNETCDF_SAFE_MODE=$j
    # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"

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

