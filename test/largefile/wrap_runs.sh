#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
outfile=`basename $1`

# disable safe mode, as large tests already run slow
export PNETCDF_SAFE_MODE=0

${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc
rm -f ${TESTOUTDIR}/$outfile.nc

# echo ""

if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
   # echo "---- testing burst buffering"
   export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
   ${TESTSEQRUN} $1              ${TESTOUTDIR}/$outfile.nc
   unset PNETCDF_HINTS

   ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc

   rm -f ${TESTOUTDIR}/$outfile.nc
fi
