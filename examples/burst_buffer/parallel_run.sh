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

# let NTHREADS=$1*6-1
NTHREADS=`expr $1 \* 6 - 1`

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for i in ${check_PROGRAMS} ; do
    # echo "---- exec=$i"
    for j in ${safe_modes} ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        # echo "${MPIRUN} ./$i -q -b ${TESTOUTDIR} ${TESTOUTDIR}/$i.nc"
        ${MPIRUN} ./$i -q -b ${TESTOUTDIR} ${TESTOUTDIR}/$i.nc
        if test $? = 0 ; then
           echo "PASS:  C  parallel run on $1 processes --------------- $i"
        fi

        # echo "--- validating file ${TESTOUTDIR}/$i.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        # echo ""
    done
done

