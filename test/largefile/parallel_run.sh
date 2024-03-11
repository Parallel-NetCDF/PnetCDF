#!/bin/sh
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# turn off safe mode for large tests
export PNETCDF_SAFE_MODE=0

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
    ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc

    # echo "--- validating file ${TESTOUTDIR}/$i.nc"
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc
    # echo ""

    if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
       # echo "---- test burst buffering feature"
       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
       ${MPIRUN} ./$i ${TESTOUTDIR}/$i.bb.nc
       unset PNETCDF_HINTS

       # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
       ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc
    fi

    rm -f ${OUTDIR}/$i.bb.nc
done

