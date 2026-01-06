#!/bin/bash
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DRY_RUN=no
VERBOSE=no

exe_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $MPIRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $MPIRUN $@
   fi
}

seq_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $TESTSEQRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $TESTSEQRUN $@
   fi
}

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do

   if test "x$i" = "xtest_inq_format" ; then
      CMD_OPTS="-q -i ${srcdir}"
   elif test "x$i" = "xtst_corrupt" ; then
      CMD_OPTS="${srcdir}"
   elif test "x$i" = "xtst_open_cdf5" ; then
      CMD_OPTS=${srcdir}/bad_begin.nc5
   else
      CMD_OPTS="-q -o ${TESTOUTDIR}/$i.nc"
   fi

   for j in ${safe_modes} ; do

      export PNETCDF_SAFE_MODE=$j
      if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
         echo "Line ${LINENO}: PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE"
      fi

      exe_cmd ./$i $CMD_OPTS

   done # safe_modes
done # check_PROGRAMS

