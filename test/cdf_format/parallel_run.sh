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
# echo "check_PROGRAMS=${check_PROGRAMS}"
# echo "srcdir = ${srcdir}"

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for j in ${safe_modes} ; do
    export PNETCDF_SAFE_MODE=$j
    # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${MPIRUN} ./test_inq_format ${srcdir}
    ${MPIRUN} ./cdf_type ${TESTOUTDIR}/cdf_type.nc
    ${MPIRUN} ./dim_cdf12 ${TESTOUTDIR}/dim_cdf12.nc

    # echo "--- validating file ${TESTOUTDIR}/cdf_type.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/cdf_type.nc
    # echo "--- validating file ${TESTOUTDIR}/dim_cdf12.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/dim_cdf12.nc

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       # echo "test burst buffering feature"
       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${MPIRUN} ./test_inq_format ${srcdir}
       ${MPIRUN} ./cdf_type ${TESTOUTDIR}/cdf_type.bb.nc
       ${MPIRUN} ./dim_cdf12 ${TESTOUTDIR}/dim_cdf12.bb.nc
       unset PNETCDF_HINTS

       # echo "--- validating file ${TESTOUTDIR}/cdf_type.bb.nc"
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/cdf_type.bb.nc
       # echo "--- validating file ${TESTOUTDIR}/dim_cdf12.bb.nc"
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/dim_cdf12.bb.nc

       # echo "--- ncmpidiff cdf_type.nc cdf_type.bb.nc ---"
       ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/cdf_type.nc ${TESTOUTDIR}/cdf_type.bb.nc
       # echo "--- ncmpidiff dim_cdf12.nc dim_cdf12.bb.nc ---"
       # ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/dim_cdf12.nc ${TESTOUTDIR}/dim_cdf12.bb.nc
    fi
done

