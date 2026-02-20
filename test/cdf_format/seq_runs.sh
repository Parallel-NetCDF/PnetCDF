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
      echo "Line $lineno CMD: $TESTSEQRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $TESTSEQRUN $@
   fi
}

exe_name=`basename $1`

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

if test "x$exe_name" = xtest_inq_format ; then
   run_cmd ./$1 -q -i ${srcdir}
elif test "x$exe_name" = xtst_open_cdf5 ; then
   # check files with corrupted header
   run_cmd ./$1 ${srcdir}/bad_begin.nc5
elif test "x$exe_name" = xtst_corrupt ; then
   # check files with corrupted header
   run_cmd ./$1 ${srcdir}
else
   run_cmd ./$1 -q -o ${TESTOUTDIR}/${exe_name}.nc
fi

