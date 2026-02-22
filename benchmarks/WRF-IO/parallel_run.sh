#!/bin/bash
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
# set -e

DRY_RUN=no
VERBOSE=no

run_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $MPIRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $MPIRUN $@
   fi
   if test $? != 0 ; then
      echo "FAIL: nprocs=$1 ---- $i $TEST_OPTS"
      exit 1
   fi
}

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

# let NTHREADS=$1*6-1
NTHREADS=`expr $1 \* 6 - 1`

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

TEST_MPIIO_MODES="0 1"

for i in ${check_PROGRAMS} ; do
    # Capture start time in seconds and nanoseconds
    start_time=$(date +%s.%1N)

    OUT_PREFIX="${TESTOUTDIR}/$i"

    for mpiio_mode in ${TEST_MPIIO_MODES} ; do
        if test "$mpiio_mode" = 1 ; then
           USEMPIO_HINTS="pnc_driver=mpiio"
           DRIVER_OUT_FILE="${OUT_PREFIX}.mpio"
           driver_hint=" MPIO"
        else
           USEMPIO_HINTS="pnc_driver=pncio"
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
        TEST_OPTS="$driver_hint $ina_hint"

        PNETCDF_HINTS=
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
        # echo "${LINENO}: PNETCDF_HINTS=$PNETCDF_HINTS"

        CMD_OPTS="-q -d -g -n 2 -y 100 -x 10 -i ${srcdir}/wrf_header.txt"
        run_cmd ./$i $CMD_OPTS -w $OUT_FILE.nc -r $OUT_FILE.nc

        # echo "${LINENO}:--- validating file ${OUT_FILE}.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "${LINENO}: ---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           # echo "${LINENO}: PNETCDF_HINTS=$PNETCDF_HINTS"
           run_cmd ./$i $CMD_OPTS -w $OUT_FILE.bb.nc -r $OUT_FILE.bb.nc

           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # echo "${LINENO}: --- validating file ${OUT_FILE}.bb.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc

           DIFF_OPT="-q"
           run_cmd ${NCMPIDIFF} $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "${LINENO}: test netCDF-4 feature"
           run_cmd ./$i $CMD_OPTS -w $OUT_FILE.nc4 -r $OUT_FILE.nc4
           # Validator does not support nc4
        fi
    done # intra_aggr
    done # mpiio_mode

    DIFF_OPT="-q"
    if test "x$TEST_MPIIO_MODES" = "x0 1" ; then
       run_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc
       run_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc
    fi
    run_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc $OUT_PREFIX.pncio.ina.nc

    rm -f ${OUTDIR}/$i*nc*

    end_time=$(date +%s.%1N)

    # Calculate difference (requires bc for floating point math)
    elapsed_time=$(echo "$end_time - $start_time" | bc)

    fixed_length=48
    printf "*** TESTING  %-${fixed_length}s   -- pass (%4ss)\n" "$i" "$elapsed_time"

done # check_PROGRAMS

