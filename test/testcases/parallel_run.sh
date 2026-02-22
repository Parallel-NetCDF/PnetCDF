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

   # # SECONDS=0
   # start_ns=$(date +%s.%4N)

   exe_name=`basename $i`

   # PNCIO driver does not support vard APIs
   if test "x$exe_name" = xtest_vardf90 || test "x$exe_name" = xtest_vardf ; then
      export PNETCDF_HINTS="pnc_driver=mpiio;$PNETCDF_HINTS"
   fi

   run_cmd ./$i -q -o ${TESTOUTDIR}/${exe_name}.nc

   # # echo "Elapsed: $SECONDS seconds"
   # end_ns=$(date +%s.%4N)
   # # Calculate difference (requires bc for floating point math)
   # elapsed_ns=$(echo "$end_ns - $start_ns" | bc)
   # echo "Elapsed time: ${elapsed_ns} seconds"

done # check_PROGRAMS

