#!/bin/bash
#
# Copyright (C) 2026, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DRY_RUN=no
VERBOSE=yes

if test "x$1" = "x./ncmpidiff" ; then
   RUN_CMD=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/4/g"`
else
   RUN_CMD=$TESTSEQRUN
fi

run_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $RUN_CMD $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $RUN_CMD $@
   fi
}

exe_name=`basename $1`

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

# ncmpidiff/ncdiff should error out when two input file names are identical.
IDENT_FILE=${srcdir}/tst_file.nc

run_cmd $1 ${IDENT_FILE} ${IDENT_FILE}

