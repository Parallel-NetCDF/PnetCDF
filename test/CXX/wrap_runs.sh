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
    export PNETCDF_SAFE_MODE=$j
    # echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${TESTSEQRUN} ./$1            ${TESTOUTDIR}/$outfile.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc
    # echo ""

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       echo ""
       echo "---- testing burst buffering"

       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.bb.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.bb.nc

       # echo "--- ncmpidiff $outfile.nc $outfile.bb.nc ---"
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$outfile.nc ${TESTOUTDIR}/$outfile.bb.nc
    fi
done
rm -f ${OUTDIR}/$outfile.nc
rm -f ${OUTDIR}/$outfile.bb.nc

