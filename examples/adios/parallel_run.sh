#!/bin/bash
#
# Copyright (C) 2019, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

let NTHREADS=$1*6-1

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    PNETCDF_HINTS=
    if test "${PNETCDF_DEBUG}" = 1 ; then # test only in safe mode
       PNETCDF_HINTS="romio_no_indep_rw=true"
    fi
    if test "x$MIMIC_LUSTRE" != x1 ; then
       PNETCDF_HINTS="cb_nodes=2;$PNETCDF_HINTS"
    fi
    export PNETCDF_HINTS=$PNETCDF_HINTS
    if test "$i" = read_metadata ; then
       ${MPIRUN} ./$i -q ${top_srcdir}/test/adios/attributes.bp
    else
       ${MPIRUN} ./$i -q ${top_srcdir}/test/adios/arrays.bp
    fi
    if test $? = 0 ; then
       echo "PASS:  C  parallel run on $1 processes --------------- $i"
    fi
done

