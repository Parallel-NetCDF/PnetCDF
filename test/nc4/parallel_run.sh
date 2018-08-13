#!/bin/sh
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

for i in ${check_PROGRAMS} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        if test "$i" = compressed ; then
           ${MPIRUN} ./$i ${srcdir}
        else
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc
        fi
    done
done

