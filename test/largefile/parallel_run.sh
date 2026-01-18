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

if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

if test "x$MIMIC_LUSTRE" != x1 ; then
   PNETCDF_HINTS="cb_nodes=2"
fi

for i in ${check_PROGRAMS} ; do

   exe_name=`basename $i`

   for j in ${safe_modes} ; do

      export PNETCDF_SAFE_MODE=$j
      if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
         echo "Line ${LINENO}: PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE"
      fi

      run_cmd ./$i -q -o ${TESTOUTDIR}/${exe_name}.nc

   done # safe_modes

done # check_PROGRAMS

