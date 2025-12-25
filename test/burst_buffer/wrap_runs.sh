#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
outfile=`basename $1`

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

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
for mpiio_mode in 0 1 ; do
    if test "$mpiio_mode" = 1 ; then
       USEMPIO_HINTS="nc_pncio=disable"
    else
       USEMPIO_HINTS="nc_pncio=enable"
    fi

for bb_mode in 1 ; do
    PNETCDF_HINTS=
    if test "x$SAFE_HINTS" != x ; then
       PNETCDF_HINTS="$SAFE_HINTS"
    fi
    if test "x$USEMPIO_HINTS" != x ; then
       PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
    fi
    if test "$bb_mode" = 1 ; then
       PNETCDF_HINTS="$PNETCDF_HINTS;nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
    fi
    export PNETCDF_HINTS="$PNETCDF_HINTS"
    export PNETCDF_SAFE_MODE=$j
    # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"

    ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
    unset PNETCDF_HINTS
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc
done
done
done

rm -f ${OUTDIR}/$outfile.nc
rm -f ${OUTDIR}/$outfile.nc_0_0.data
rm -f ${OUTDIR}/$outfile.nc_0_0.meta

