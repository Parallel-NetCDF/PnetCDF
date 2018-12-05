#!/bin/sh
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# let NTHREADS=$1*6-1
NTHREADS=`expr $1 \* 6 - 1`

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        if test $i = "pnetcdf-read-from-master" ; then
           ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-from-master.nc
        elif test $i = "pnetcdf-read-nfiles" ; then
           ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-nfiles.nc
        elif test $i = "pnetcdf-read-standard" ; then
           ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-standard.nc
        elif test $i = "pnetcdf-read-flexible" ; then
           ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-flexible.nc
        elif test $i = "pnetcdf-read-nb" ; then
           ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-nb.nc
        else
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc
        fi
        if test $? = 0 ; then
           if test $i = "pnetcdf-write-bufferedf77" ; then
              echo "PASS: F77 parallel run on $1 processes --------------- $i"
           elif test $i = "pnetcdf-write-bufferedf" ; then
              echo "PASS: F90 parallel run on $1 processes --------------- $i"
           else
              echo "PASS:  C  parallel run on $1 processes --------------- $i"
           fi
        fi

        if test "$i" = pthread ; then
           # each MPI process created 6 threads
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.$k
           done
           continue
        fi

        if test $i = "pnetcdf-write-nfiles" ; then
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.0-4.nc
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.1-4.nc
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.2-4.nc
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.3-4.nc
        elif test $i != "pnetcdf-read-from-master" -a \
                  $i != "pnetcdf-read-nfiles" -a \
                  $i != "pnetcdf-read-standard" -a \
                  $i != "pnetcdf-read-flexible" -a \
                  $i != "pnetcdf-read-nb" ; then
           # echo "--- validating file ${TESTOUTDIR}/$i.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        fi
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "test burst buffering feature"
           export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           if test $i = "pnetcdf-read-from-master" ; then
              ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-from-master.bb.nc
           elif test $i = "pnetcdf-read-nfiles" ; then
              ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-nfiles.bb.nc
           elif test $i = "pnetcdf-read-standard" ; then
              ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-standard.bb.nc
           elif test $i = "pnetcdf-read-flexible" ; then
              ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-flexible.bb.nc
           elif test $i = "pnetcdf-read-nb" ; then
              ${MPIRUN} ./$i ${TESTOUTDIR}/pnetcdf-write-nb.bb.nc
           else
              ${MPIRUN} ./$i ${TESTOUTDIR}/$i.bb.nc
           fi
           if test $? = 0 ; then
              if test $i = "pnetcdf-write-bufferedf77" ; then
                 echo "PASS: F77 parallel run on $1 processes --------------- $i"
              elif test $i = "pnetcdf-write-bufferedf" ; then
                 echo "PASS: F90 parallel run on $1 processes --------------- $i"
              else
                 echo "PASS:  C  parallel run on $1 processes --------------- $i"
              fi
           fi
           unset PNETCDF_HINTS

           if test $i = "pnetcdf-write-nfiles" ; then
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.0-4.nc
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.1-4.nc
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.2-4.nc
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc.3-4.nc
           elif test $i != "pnetcdf-read-from-master" -a \
                     $i != "pnetcdf-read-nfiles" -a \
                     $i != "pnetcdf-read-standard" -a \
                     $i != "pnetcdf-read-flexible" -a \
                     $i != "pnetcdf-read-nb" ; then
              # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc

              # echo "--- ncmpidiff $i.nc $i.bb.nc ---"
              ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc ${TESTOUTDIR}/$i.bb.nc
           fi
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc4 4
           # Validator does not support nc4
        fi
    done
done

