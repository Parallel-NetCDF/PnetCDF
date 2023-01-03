#!/bin/sh
#
# Copyright (C) 2022, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

NPROCS=4
if test "x$1" != x ; then
   NPROCS=$1
fi
OUTDIR=TESTOUTDIR
MPIRUN="TESTMPIRUN"
MPIRUN=`echo ${MPIRUN} | SED_CMD -e "s/NP/${NPROCS}/g"`
run_BURST_BUFFER=ENABLE_BURST_BUFFER
run_NETCDF4=ENABLE_NETCDF4

for i in check_PROGRAMS ; do
    if test $i = get_vara ; then
       # get_vara reads the file 'put_vara.nc' created by put_vara
       ${MPIRUN} ./$i -q ${OUTDIR}/put_vara.nc
    else
       ${MPIRUN} ./$i -q ${OUTDIR}/$i.nc
    fi
    if test $? = 0 ; then
       echo "PASS: F77 parallel run on ${NPROCS} processes --------------- $i"
    fi

    if test "x${run_BURST_BUFFER}" = x1 ; then
       # echo "test burst buffering feature"
       export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${OUTDIR};nc_burst_buf_overwrite=enable"
       if test $i = get_vara ; then
          ${MPIRUN} ./$i -q ${OUTDIR}/put_vara.bb.nc
       else
          ${MPIRUN} ./$i -q ${OUTDIR}/$i.bb.nc
       fi
       if test $? = 0 ; then
          echo "PASS: F77 parallel run on ${NPROCS} processes --------------- $i"
       fi
       unset PNETCDF_HINTS
    fi

    if test "x${run_NETCDF4}" = x1 ; then
       # echo "test netCDF-4 feature"
       ${MPIRUN} ./$i ${OUTDIR}/$i.nc4 4
    fi
    # delete output file
    rm -f ${OUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.bb.nc
done

