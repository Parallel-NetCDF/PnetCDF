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

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${TESTPROGRAMS} ; do
    for j in ${safe_modes} ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./$i -f ${TESTOUTDIR}/$i.nc -s 2

        # echo "--- validating file ${TESTOUTDIR}/$i.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        # echo "--- validating file ${TESTOUTDIR}/$i.nc.subfile_0.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.subfile_0.nc
        # echo ""

        # skip burst buffering test, as it has not supported subfiling yet
        continue

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           echo ""
           echo "---- testing burst buffering"

           export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${TESTSEQRUN} ./$i -f ${TESTOUTDIR}/$i.bb.nc -s 2
           unset PNETCDF_HINTS

           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc
           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc

           # echo "--- ncmpidiff $i.nc.subfile_0.nc $i.bb.nc.subfile_0.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc.subfile_0.nc ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc
        fi
    done
    rm -f ${OUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc.subfile_0.nc
    rm -f ${OUTDIR}/$i.bb.nc
    rm -f ${OUTDIR}/$i.bb.nc.subfile_0.nc
done

