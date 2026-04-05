#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

outfile=`basename $1`

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for mpiio_mode in 0 1 ; do
    if test "$mpiio_mode" = 1 ; then
       USEMPIO_HINTS="nc_driver=mpiio"
    else
       USEMPIO_HINTS="nc_driver=pncio"
    fi

    PNETCDF_HINTS=
    if test "x$USEMPIO_HINTS" != x ; then
       PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
    fi

    export PNETCDF_HINTS="$PNETCDF_HINTS"
    # echo "PNETCDF_HINTS=$PNETCDF_HINTS"

    ${TESTSEQRUN} $1 ${TESTOUTDIR}/$outfile.nc
done

rm -f ${OUTDIR}/$outfile.nc
rm -f ${OUTDIR}/$outfile.nc.cdf4

