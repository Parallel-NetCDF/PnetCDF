#!/bin/sh
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
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
    for intra_aggr in 0 1 ; do
        if test "$j" = 1 ; then # test only in safe mode
           export PNETCDF_HINTS="romio_no_indep_rw=true"
        else
           export PNETCDF_HINTS=
        fi
        if test "$intra_aggr" = 1 ; then
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_num_aggrs_per_node=2"
        fi
        # echo "PNETCDF_HINTS=${PNETCDF_HINTS}"

        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        OPTS="-y 100 -x 100 -i ${srcdir}/wrf_header.txt"
        OPTS="$OPTS -w ${TESTOUTDIR}/$i.nc -r ${TESTOUTDIR}/$i.nc"
        # echo "${MPIRUN} ./$i -q ${OPTS}"
        ${MPIRUN} ./$i -q ${OPTS}
        if test $? = 0 ; then
           echo "PASS:  C  parallel run on $1 processes --------------- $i"
        fi

        unset PNETCDF_HINTS
        # echo "--- validating file ${TESTOUTDIR}/$i.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        # echo ""
    done
    done
    rm -f ${OUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.nc4
done

