#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

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

    ${TESTSEQRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc
    ${TESTSEQRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc
    # echo "--- validating file ${TESTOUTDIR}/pres_temp_4D.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q   ${TESTOUTDIR}/pres_temp_4D.nc
    # echo ""

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       echo ""
       echo "---- testing burst buffering"
       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${TESTSEQRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.bb.nc
       ${TESTSEQRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.bb.nc
       unset PNETCDF_HINTS
       # echo "--- validating file ${TESTOUTDIR}/pres_temp_4D.bb.nc"
       ${TESTSEQRUN} ${VALIDATOR} -q   ${TESTOUTDIR}/pres_temp_4D.bb.nc
       # echo ""

       # echo "--- ncmpidiff pres_temp_4D.nc pres_temp_4D.bb.nc ---"
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/pres_temp_4D.nc ${TESTOUTDIR}/pres_temp_4D.bb.nc
    fi
    # echo ""

    if test "x${ENABLE_NETCDF4}" = x1 ; then
        ${TESTSEQRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc4 4
        ${TESTSEQRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc4
        # Validator does not support nc4
    fi
    # echo ""
done
done
rm -f ${OUTDIR}/*.nc
rm -f ${OUTDIR}/*.nc4
