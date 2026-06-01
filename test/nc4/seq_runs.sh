#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

outfile=`basename $1`

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

if test "x$ENABLE_GIO" = x0 ; then
   IO_MODES="mpiio"
else
   IO_MODES="gio mpiio"
   if test "x$GIO_ONLY" = x1 ; then
      IO_MODES="gio"
   fi
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

if test "x$GIO_ONLY" = x1 ; then
   export PNETCDF_HINTS="nc_driver=gio"
fi

for io_mode in $IO_MODES ; do
    if test "x$io_mode" = xmpiio ; then
       USEMPIO_HINTS="nc_driver=mpiio"
    else
       USEMPIO_HINTS="nc_driver=gio"
    fi

    PNETCDF_HINTS=
    if test "x$USEMPIO_HINTS" != x ; then
       PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
    fi

    export PNETCDF_HINTS="$PNETCDF_HINTS"
    # echo "PNETCDF_HINTS=$PNETCDF_HINTS"

    ${TESTSEQRUN} $1 ${TESTOUTDIR}/$outfile.nc
done

rm -f ${OUTDIR}/$outfile.nc
rm -f ${OUTDIR}/$outfile.nc.cdf4

