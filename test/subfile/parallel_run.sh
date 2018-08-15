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

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${MPIRUN} ./$i -f ${TESTOUTDIR}/$i.nc -s 2

        # echo "--- validating file ${TESTOUTDIR}/$i.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        # echo "--- validating file ${TESTOUTDIR}/$i.nc.subfile_0.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.subfile_0.nc
        # echo "--- validating file ${TESTOUTDIR}/$i.nc.subfile_1.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.subfile_1.nc
        # echo ""

        # skip burst buffering test, as it has not supported subfiling yet
        continue

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "test burst buffering feature"
           export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${MPIRUN} ./$i -f ${TESTOUTDIR}/$i.bb.nc -s 2
           unset PNETCDF_HINTS

           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc
           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc
           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc.subfile_1.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.subfile_1.nc

           # echo "--- ncmpidiff $i.nc $i.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc ${TESTOUTDIR}/$i.bb.nc
           # echo "--- ncmpidiff $i.nc.subfile_0.nc $i.bb.nc.subfile_0.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc.subfile_0.nc ${TESTOUTDIR}/$i.bb.nc.subfile_0.nc
           # echo "--- ncmpidiff $i.nc.subfile_1.nc $i.bb.nc.subfile_1.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc.subfile_1.nc ${TESTOUTDIR}/$i.bb.nc.subfile_1.nc
        fi
    done
done

