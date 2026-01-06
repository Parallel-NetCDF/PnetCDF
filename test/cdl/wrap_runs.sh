#!/bin/bash
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DRY_RUN=no
VERBOSE=no

seq_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $TESTSEQRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $TESTSEQRUN $@
   fi
}

outfile=`basename $1`

CMD_OPTS=
if test "x$1" = "x./tst_cdl_hdr_parser" ; then
   CMD_OPTS="-q -o ${TESTOUTDIR}/$outfile.nc -i ${srcdir}/cdl_header.txt"
fi

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for j in ${safe_modes} ; do

   # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE"
   PNETCDF_HINTS="$SAFE_HINTS"

   seq_cmd $1 $CMD_OPTS

done
