#!/bin/bash
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

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

    # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE"
    export PNETCDF_SAFE_MODE=$j

    ${TESTSEQRUN} ./pres_temp_4D_wr_rd -q -o ${TESTOUTDIR}/pres_temp_4D.nc
    # echo ""

done
