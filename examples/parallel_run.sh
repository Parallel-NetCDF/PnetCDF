#!/bin/bash
#
# Copyright (C) 2018, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DRY_RUN=no
VERBOSE=no

exe_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}

   cmd=`basename $1`
   if test "x$MIMIC_LUSTRE" = x1 && test "x$cmd" = xncmpidiff ; then
      # echo "export MIMIC_STRIPE_SIZE=1048576"
      export MIMIC_STRIPE_SIZE=1048576
   fi

   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $MPIRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $MPIRUN $@
   fi

   if test "x$MIMIC_LUSTRE" = x1 && test "x$cmd" = xncmpidiff ; then
      # echo "unset MIMIC_STRIPE_SIZE"
      unset MIMIC_STRIPE_SIZE
   fi
}

seq_cmd() {
   local lineno=${BASH_LINENO[$((${#BASH_LINENO[@]} - 2))]}
   if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
      echo "Line $lineno CMD: $TESTSEQRUN $@"
   fi
   if test "x$DRY_RUN" = xno ; then
      $TESTSEQRUN $@
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

# echo "${LINENO}: PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
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
           safe_hint="  SAFE"
        else
           safe_hint="NOSAFE"
        fi
        OUT_PREFIX="${TESTOUTDIR}/$i"

    for no_indep_rw in true false ; do
        no_indep_rw_hint="romio_no_indep_rw=$no_indep_rw"

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

        PNETCDF_HINTS=$no_indep_rw_hint
        if test "x$SAFE_HINTS" != x ; then
           PNETCDF_HINTS="$SAFE_HINTS;$PNETCDF_HINTS"
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
        if test "x$VERBOSE" = xyes || test "x$DRY_RUN" = xyes ; then
           echo "Line ${LINENO}: PNETCDF_SAFE_MODE=$PNETCDF_SAFE_MODE PNETCDF_HINTS=$PNETCDF_HINTS"
        fi

        TEST_OPTS="$safe_hint $driver_hint $ina_hint"

        CMD_OPT=-q
        IN_FILE=
        if test "$i" = create_from_cdl ; then
           IN_FILE=${srcdir}/cdl_header.txt
        fi

        if test "$i" = pthread ; then
           # each MPI process created 6 threads
           exe_cmd ./$i $CMD_OPT ${OUT_FILE}.nc
           for k in `seq 0 ${NTHREADS}` ; do
               seq_cmd ${VALIDATOR} -q ${OUT_FILE}.nc.$k
           done
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
           continue
        elif test "$i" = put_vara ; then
           exe_cmd ./$i $CMD_OPT ${OUT_FILE}.nc
           seq_cmd ${VALIDATOR} -q ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi

           exe_cmd ./get_vara $CMD_OPT ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- get_vara\n" $2 "$TEST_OPTS"
           fi
        elif test "$i" = get_vara ; then
           continue
        elif test "$i" = create_from_cdl ; then
           # create_from_cdl reads a CDL header file
           exe_cmd ./$i $CMD_OPT -o ${OUT_FILE}.nc $IN_FILE
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
        else
           exe_cmd ./$i $CMD_OPT ${OUT_FILE}.nc
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS"
           fi
        fi

        seq_cmd ${VALIDATOR} -q ${OUT_FILE}.nc

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           if test "$i" = create_from_cdl ; then
              # create_from_cdl reads a CDL header file
              exe_cmd ./$i -q -o ${OUT_FILE}.bb.nc $IN_FILE
           else
              exe_cmd ./$i $CMD_OPT ${OUT_FILE}.bb.nc
           fi
           if test $? = 0 ; then
              printf "PASS: %-3s nprocs=$1 %-${fixed_length}s   -------- $i\n" $2 "$TEST_OPTS BB"
           fi
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           seq_cmd ${VALIDATOR} -q ${OUT_FILE}.bb.nc

           # compare file header only for large file tests
           DIFF_OPT="-q"
           if test "$i" = create_from_cdl ; then
              DIFF_OPT+=" -h"
           fi
           exe_cmd ${NCMPIDIFF} $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           exe_cmd ./$i ${OUT_FILE}.nc4 4
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
          exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j
          exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j
          exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j
       done
    else
       exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.mpio.ina.nc
       exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc $OUT_PREFIX.pncio.nc
       exe_cmd $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc $OUT_PREFIX.pncio.ina.nc
    fi

    done # no_indep_rw
    done # safe_modes
    rm -f ${OUTDIR}/$i*nc*
done # check_PROGRAMS

