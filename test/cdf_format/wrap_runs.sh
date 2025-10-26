#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

outfile=`basename $1`

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
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

    PNETCDF_HINTS=
    if test "x$SAFE_HINTS" != x ; then
       PNETCDF_HINTS="$SAFE_HINTS"
    fi
    if test "x$USEMPIO_HINTS" != x ; then
       PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
    fi

    export PNETCDF_HINTS="$PNETCDF_HINTS"
    export PNETCDF_SAFE_MODE=$j
    # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"

    ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       echo ""
       echo "---- testing burst buffering"

       export PNETCDF_HINTS="$PNETCDF_HINTS;nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.bb.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.bb.nc

       # running ncmpidiff on dim_cdf12 on one process requires more than 2 GB
       # memory. Use more processes to check. Disable the check for now.
       if test "$1" != "./dim_cdf12" ; then
          # echo "--- ncmpidiff $outfile.nc $outfile.bb.nc ---"
          ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$outfile.nc ${TESTOUTDIR}/$outfile.bb.nc
       fi
    fi
done
done
rm -f ${OUTDIR}/$outfile.nc
rm -f ${OUTDIR}/$outfile.bb.nc
