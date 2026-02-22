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

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

if test "${PNETCDF_DEBUG}" = 1 ; then # test only in safe mode
   SAFE_HINTS="romio_no_indep_rw=true"
else
   SAFE_HINTS="romio_no_indep_rw=false"
fi

for mpiio_mode in 0 1 ; do
    if test "$mpiio_mode" = 1 ; then
       USEMPIO_HINTS="pnc_driver=mpiio"
    else
       USEMPIO_HINTS="pnc_driver=pncio"
    fi

    PNETCDF_HINTS=
    if test "x$SAFE_HINTS" != x ; then
       PNETCDF_HINTS="$SAFE_HINTS"
    fi
    if test "x$USEMPIO_HINTS" != x ; then
       PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
    fi

    export PNETCDF_HINTS="$PNETCDF_HINTS"
    # echo "PNETCDF_HINTS=$PNETCDF_HINTS"

    ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc
    # echo ""

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       echo ""
       echo "---- testing burst buffering"

       export PNETCDF_HINTS="$PNETCDF_HINTS;nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.bb.nc
       unset PNETCDF_HINTS
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.bb.nc

       # echo "--- ncmpidiff $outfile.nc $outfile.bb.nc ---"
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$outfile.nc ${TESTOUTDIR}/$outfile.bb.nc
    fi
    rm -f ${OUTDIR}/$outfile.nc
    rm -f ${OUTDIR}/$outfile.bb.nc
    rm -f ${OUTDIR}/$outfile.nc.2
    rm -f ${OUTDIR}/$outfile.bb.nc.2
done

