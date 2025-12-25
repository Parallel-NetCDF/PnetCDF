#!/bin/bash
#
# Copyright (C) 2025, Northwestern University and Argonne National Laboratory
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

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do

   CMD_OPTS=${TESTOUTDIR}/$i.nc
   if test $i = "tst_cdl_hdr_parser" ; then
      CMD_OPTS="-q -o ${TESTOUTDIR}/$i.nc ${srcdir}/cdl_header.txt"
   fi

    for j in ${safe_modes} ; do
        if test "$j" = 1 ; then # test only in safe mode
           export PNETCDF_HINTS="romio_no_indep_rw=true"
        else
           export PNETCDF_HINTS=
        fi
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${MPIRUN} ./$i ${CMD_OPTS}

    done
    rm -f ${OUTDIR}/$i.nc
done

