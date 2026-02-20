#!/bin/bash
#
# Copyright (C) 2026, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DRY_RUN=no
VERBOSE=no

run_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $MPIRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $MPIRUN $@
   fi
}

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

if test "x$MIMIC_LUSTRE" != x1 ; then
   PNETCDF_HINTS="cb_nodes=2"
fi

for i in ${check_PROGRAMS} ; do

   exe_name=`basename $i`

   if test "x$exe_name" = "xtest_inq_format" ; then
      run_cmd ./$i -q -i ${srcdir}
   elif test "x$exe_name" = "xtst_corrupt" ; then
      run_cmd ./$i ${srcdir}
   elif test "x$exe_name" = "xtst_open_cdf5" ; then
      run_cmd ./$i ${srcdir}/bad_begin.nc5
   else
      run_cmd ./$i -q -o ${TESTOUTDIR}/${exe_name}.nc
   fi

done # check_PROGRAMS

