#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for j in ${safe_modes} ; do
    export PNETCDF_SAFE_MODE=$j
    # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${TESTSEQRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc
    ${TESTSEQRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc
    # echo "--- validating file ${TESTOUTDIR}/pres_temp_4D.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q   ${TESTOUTDIR}/pres_temp_4D.nc
    # echo ""

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       # echo "test burst buffering feature"
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

    if test "${ENABLE_NETCDF4}" = 1 ; then
        ${TESTSEQRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc4 4
        ${TESTSEQRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc4
        # Validator does not support nc4
    fi
    # echo ""
done

