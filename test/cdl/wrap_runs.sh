#!/bin/sh
#
# Copyright (C) 2025, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

outfile=`basename $1`

CMD_OPTS=
if test $outfile = "tst_cdl_hdr_parser" ; then
   CMD_OPTS="-q -o ${TESTOUTDIR}/$outfile.nc ${srcdir}/cdl_header.txt"
fi

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
    ${TESTSEQRUN} $1 $CMD_OPTS
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$outfile.nc

done
rm -f ${OUTDIR}/$outfile.nc

