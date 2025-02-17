#!/bin/bash
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

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

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
        OUT_PREFIX="${TESTOUTDIR}/$i"

    for mpiio_mode in 0 1 ; do
        if test "$mpiio_mode" = 1 ; then
           USEMPIO_HINTS="nc_pncio=disable"
           DRIVER_OUT_FILE="${OUT_PREFIX}.mpio"
        else
           USEMPIO_HINTS="nc_pncio=enable"
           DRIVER_OUT_FILE="${OUT_PREFIX}.pncio"
        fi
    for intra_aggr in 0 1 ; do
        if test "$intra_aggr" = 1 ; then
           INA_HINTS="nc_num_aggrs_per_node=2"
           INA_OUT_FILE="${DRIVER_OUT_FILE}.ina"
        else
           INA_HINTS="nc_num_aggrs_per_node=0"
           INA_OUT_FILE="${DRIVER_OUT_FILE}"
        fi

        OUT_FILE=$INA_OUT_FILE

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

        if test "$i" = tst_version ; then
           ${MPIRUN} ./tst_version
           continue
        fi

        if test "$i" = tst_pthread ; then
           # each MPI process created 6 threads
           ${MPIRUN} ./tst_pthread ${OUT_FILE}.nc
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc.$k
               rm -f ${OUTDIR}/tst_pthread.nc.$k
           done
           continue
        fi

        ${MPIRUN} ./$i ${OUT_FILE}.nc

        # put_all_kinds and iput_all_kinds output 3 files
        if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
           for k in 1 2 5 ; do
               # echo "--- validating file ${OUT_FILE}.nc$k"
               ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc$k
           done
        else
           # echo "--- validating file ${OUT_FILE}.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc
        fi
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${MPIRUN} ./$i ${OUT_FILE}.bb.nc
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # put_all_kinds and iput_all_kinds output 3 files
           if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
              for k in 1 2 5 ; do
                  # echo "--- validating file ${OUT_FILE}.bb.nc$k"
                  ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc$k
                  # echo "--- ncmpidiff $i.nc$k $i.bb.nc$k ---"
                  ${MPIRUN} ${NCMPIDIFF} -q ${OUT_FILE}.nc$k ${OUT_FILE}.bb.nc$k
              done
              continue
           else
              # echo "--- validating file ${OUT_FILE}.bb.nc"
              ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc
           fi

           # skip ncmpidiff for large file
           if test "$i" = last_large_var ; then
              continue
           fi

           # echo "--- ncmpidiff $i.nc $i.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${OUT_FILE}.nc ${OUT_FILE}.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${OUT_FILE}.nc4 4
           # Validator does not support nc4
        fi
    done # intra_aggr
    done # mpiio_mode

    if test "$i" = tst_version ; then
       # this program creates no output file
       continue
    fi
    if [[ "$i" == *"vard"* ]] ; then
       continue
    fi

    DIFF_OPT="-q"
    if test "$i" = last_large_var ; then
       DIFF_OPT+=" -h"
    fi
    if test "$i" = put_all_kinds || test "$i" = iput_all_kinds ; then
       for j in 1 2 5; do
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.mpio.ina.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.mpio.ina.nc$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.nc$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.ina.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.ina.nc$j
       done
    elif test "$i" = tst_pthread ; then
       for j in `seq 0 ${NTHREADS}` ; do
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j
       done
    else
       # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc
       # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc
       # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.ina.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.ina.nc
    fi

    done # safe_modes
    rm -f ${OUTDIR}/$i*nc*
done # check_PROGRAMS

