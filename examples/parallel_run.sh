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

fixed_length=23

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
        if test "$j" = 1 ; then # test only in safe mode
           SAFE_HINTS="romio_no_indep_rw=true"
           safe_hint="  SAFE"
        else
           SAFE_HINTS="romio_no_indep_rw=false"
           safe_hint="NOSAFE"
        fi
        OUT_PREFIX="${TESTOUTDIR}/$i"

    for mpiio_mode in 0 1 ; do
        if test "$mpiio_mode" = 1 ; then
           USEMPIO_HINTS="nc_pncio=disable"
           DRIVER_OUT_FILE="${OUT_PREFIX}.mpio"
           driver_hint=" MPIO"
        else
           USEMPIO_HINTS="nc_pncio=enable"
           DRIVER_OUT_FILE="${OUT_PREFIX}.pncio"
           driver_hint="PNCIO"
        fi
    for intra_aggr in 0 1 ; do
        if test "$intra_aggr" = 1 ; then
           INA_HINTS="nc_num_aggrs_per_node=2"
           INA_OUT_FILE="${DRIVER_OUT_FILE}.ina"
           ina_hint="  INA"
        else
           INA_HINTS="nc_num_aggrs_per_node=0"
           INA_OUT_FILE="${DRIVER_OUT_FILE}"
           ina_hint="NOINA"
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
        if test "x$MIMIC_LUSTRE" != x1 ; then
           PNETCDF_HINTS="cb_nodes=2;$PNETCDF_HINTS"
        fi

        export PNETCDF_HINTS="$PNETCDF_HINTS"
        export PNETCDF_SAFE_MODE=$j
        # echo "PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"

        TEST_OPTS="$safe_hint $driver_hint $ina_hint"

        CMD_OPT=-q
        IN_FILE=
        if test "$i" = create_from_cdl ; then
           IN_FILE=${srcdir}/cdl_header.txt
        fi

        if test "$i" = pthread ; then
           # each MPI process created 6 threads
           ${MPIRUN} ./$i $CMD_OPT ${OUT_FILE}.nc
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc.$k
           done
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
           continue
        elif test "$i" = put_vara ; then
           ${MPIRUN} ./$i $CMD_OPT ${OUT_FILE}.nc
           # echo "--- validating file ${OUT_FILE}.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi

           ${MPIRUN} ./get_vara $CMD_OPT ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- get_vara\n" $2 "$TEST_OPTS"
           fi
        elif test "$i" = get_vara ; then
           continue
        elif test "$i" = create_from_cdl ; then
           # create_from_cdl reads a CDL header file
           ${MPIRUN} ./$i $CMD_OPT -o ${OUT_FILE}.nc $IN_FILE
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
        else
           # echo "${LINENO}: ${MPIRUN} ./$i $CMD_OPT ${OUT_FILE}.nc"
           ${MPIRUN} ./$i $CMD_OPT ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
        fi

        # echo "${LINENO}:--- validating file ${OUT_FILE}.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           # echo "${LINENO}:--- ./$i $CMD_OPT ${OUT_FILE}.bb.nc"
           if test "$i" = create_from_cdl ; then
              # create_from_cdl reads a CDL header file
              ${MPIRUN} ./$i -q -o ${OUT_FILE}.bb.nc $IN_FILE
           else
              ${MPIRUN} ./$i $CMD_OPT ${OUT_FILE}.bb.nc
           fi
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS BB"
           fi
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # echo "--- validating file ${OUT_FILE}.bb.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc

           # compare file header only for large file tests
           DIFF_OPT="-q"
           if test "$i" = create_from_cdl ; then
              DIFF_OPT+=" -h"
           fi
           # echo "${LINENO}: --- ncmpidiff $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${OUT_FILE}.nc4 4
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
           # Validator does not support nc4
        fi
    done # intra_aggr
    done # mpiio_mode

    if [[ "$i" == *"vard"* ]] ; then
       continue
    fi

    if test "$i" = get_vara ; then
       continue
    fi

    DIFF_OPT="-q"
    if test "$i" = create_from_cdl ; then
       DIFF_OPT+=" -h"
    fi
    if test "$i" = pthread ; then
       for j in `seq 0 ${NTHREADS}` ; do
          # echo "${LINENO}: --- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.pncio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j
       done
    else
       # echo "${LINENO}: --- ncmpidiff $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc
       # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc
       # echo "--- ncmpidiff $OUT_PREFIX.pncio.nc $OUT_PREFIX.pncio.ina.nc ---"
       $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc $OUT_PREFIX.pncio.ina.nc
    fi

    done # safe_modes
    rm -f ${OUTDIR}/$i*nc*
done # check_PROGRAMS

