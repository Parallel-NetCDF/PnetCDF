#!/bin/sh
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"

for j in 0 1 ; do
    export PNETCDF_SAFE_MODE=$j
    # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${MPIRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc
    ${MPIRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc
    # echo "--- validating file ${TESTOUTDIR}/pres_temp_4D.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/pres_temp_4D.nc
    # echo ""
done

if test "x${TESTBB}" = x1 ; then
   # echo "test burst buffering feature"
   export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
   ${MPIRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.bb.nc
   ${MPIRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.bb.nc
   unset PNETCDF_HINTS

   # echo "--- validating file ${TESTOUTDIR}/pres_temp_4D.bb.nc"
   ${TESTSEQRUN} ${VALIDATOR} -q   ${TESTOUTDIR}/pres_temp_4D.bb.nc

   # echo "--- ncmpidiff pres_temp_4D.nc pres_temp_4D.bb.nc ---"
   ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/pres_temp_4D.nc ${TESTOUTDIR}/pres_temp_4D.bb.nc
fi

if test "x${TESTNETCDF4}" = x1 ; then
   # echo "test netCDF-4 feature"
   ${MPIRUN} ./pres_temp_4D_wr ${TESTOUTDIR}/pres_temp_4D.nc4 4
   ${MPIRUN} ./pres_temp_4D_rd ${TESTOUTDIR}/pres_temp_4D.nc4
   # Validator does not support nc4
fi

