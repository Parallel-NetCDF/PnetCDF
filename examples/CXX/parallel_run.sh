#!/bin/bash
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# let NTHREADS=$1*6-1
NTHREADS=`expr $1 \* 6 - 1`

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
        if test "$j" = 1 ; then # test only in safe mode
           SAFE_HINTS="romio_no_indep_rw=true"
        else
           SAFE_HINTS="romio_no_indep_rw=false"
        fi
    for mpiio_mode in 0 1 ; do
        if test "$mpiio_mode" = 1 ; then
           USEMPIO_HINTS="nc_pncio=disable"
        else
           USEMPIO_HINTS="nc_pncio=enable"
        fi
    for intra_aggr in 0 1 ; do
        if test "$intra_aggr" = 1 ; then
           INA_HINTS="nc_num_aggrs_per_node=2"
        else
           INA_HINTS="nc_num_aggrs_per_node=0"
        fi

        if [[ "$i" == *"vard"* ]] ; then
           if test "x$mpiio_mode" == x0 || test "x$intra_aggr" == x1 ; then
              # vard APIs are not supported when using PNCIO
              continue
           fi
        fi

        PNETCDF_HINTS=
        if test "x$SAFE_HINTS" != x ; then
           PNETCDF_HINTS="$SAFE_HINTS"
        fi
        if test "x$USEMPIO_HINTS" != x ; then
           PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
        fi
        if test "x$INA_HINTS" != x ; then
           PNETCDF_HINTS="$INA_HINTS;$PNETCDF_HINTS"
        fi

        export PNETCDF_HINTS="$PNETCDF_HINTS"
        export PNETCDF_SAFE_MODE=$j
	    # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"

        if test $i = get_vara ; then
           ${MPIRUN} ./$i -q ${TESTOUTDIR}/put_vara.nc
        else
           ${MPIRUN} ./$i -q ${TESTOUTDIR}/$i.nc
        fi
        if test $? = 0 ; then
           echo "PASS: CXX parallel run on $1 processes --------------- $i"
        fi

        if test "$i" = pthread ; then
           # each MPI process created 6 threads
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc.$k
           done
           continue
        fi

        if test $i != get_vara ; then
           # echo "--- validating file ${TESTOUTDIR}/$i.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        fi
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           if test $i = get_vara ; then
              ${MPIRUN} ./$i -q ${TESTOUTDIR}/put_vara.bb.nc
           else
              ${MPIRUN} ./$i -q ${TESTOUTDIR}/$i.bb.nc
           fi
           if test $? = 0 ; then
              echo "PASS: CXX parallel run on $1 processes --------------- $i"
           fi
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           if test $i != get_vara ; then
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
    done
    # delete output file
    if test $i = get_vara ; then
       rm -f ${OUTDIR}/put_vara.nc
       rm -f ${OUTDIR}/put_vara.bb.nc
    elif test $i != put_vara ; then
       rm -f ${OUTDIR}/$i.nc
       rm -f ${OUTDIR}/$i.bb.nc
    fi
done

