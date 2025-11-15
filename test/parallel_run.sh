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
    if test "$i" = tst_io ; then
       # this is designed to run 1 process
       continue
    fi
    if test "$i" = tst_version ; then
       # this program read only and creates no output file
       ${MPIRUN} ./$i
       continue
    fi
    if test "$i" = tst_open_cdf5 ; then
       # this program read only and creates no output file
       ${MPIRUN} ./$i ${srcdir}/bad_begin.nc5
       continue
    fi
    if test "$i" = tst_corrupt ; then
       # this program read only and creates no output file
       ${MPIRUN} ./$i ${srcdir}
       continue
    fi

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
        TEST_OPTS="$safe_hint $driver_hint $ina_hint"

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

        if test "$i" = tst_pthread ; then
           # each MPI process created 6 threads
           ${MPIRUN} ./$i ${OUT_FILE}.nc
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc.$k
               rm -f ${OUTDIR}/tst_pthread.nc.$k
           done
           continue
        elif test "$i" = pres_temp_4D_wr ; then
           ${MPIRUN} ./$i ${OUT_FILE}.nc
           # echo "--- validating file ${OUT_FILE}.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc

           ${MPIRUN} ./pres_temp_4D_rd ${OUT_FILE}.nc
        elif test "$i" = pres_temp_4D_rd ; then
           continue
        elif test "$i" = test_inq_format ; then
           ${MPIRUN} ./$i ${srcdir}
           continue
        elif test "$i" = "tst_cdl_hdr_parser" ; then
           ${MPIRUN} ./$i -q -o ${OUT_FILE}.nc ${srcdir}/cdl_header.txt
           continue
        elif test "$i" = mcoll_perf ; then
           # echo "${LINENO}: ${MPIRUN} ./$i ${OUT_FILE}"
           ${MPIRUN} ./$i ${OUT_FILE}
        else
           # echo "${LINENO}: ${MPIRUN} ./$i ${OUT_FILE}.nc"
           ${MPIRUN} ./$i ${OUT_FILE}.nc
        fi

        # put_all_kinds and iput_all_kinds output 3 files
        if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
           for k in 1 2 5 ; do
               # echo "--- validating file ${OUT_FILE}.nc$k"
               ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc$k
           done
        elif test "$i" = mcoll_perf ; then
           for j in `seq 0 9` ; do
              ext="2.4.$j.nc"
              # echo "${LINENO}:--- validating file ${OUT_FILE}.$ext"
              ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.$ext
           done
        else
           # echo "${LINENO}:--- validating file ${OUT_FILE}.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.nc
        fi

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           if test "$i" = mcoll_perf ; then
              # echo "${LINENO}:--- ${MPIRUN} ./$i ${OUT_FILE}.bb"
              ${MPIRUN} ./$i ${OUT_FILE}.bb
           else
              ${MPIRUN} ./$i ${OUT_FILE}.bb.nc
           fi
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # put_all_kinds and iput_all_kinds output 3 files
           if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
              for k in 1 2 5 ; do
                  # echo "--- validating file ${OUT_FILE}.bb.nc$k"
                  ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc$k
                  # echo "--- ncmpidiff ${OUT_FILE}.nc$k ${OUT_FILE}.bb.nc$k ---"
                  ${MPIRUN} ${NCMPIDIFF} -q ${OUT_FILE}.nc$k ${OUT_FILE}.bb.nc$k
              done
              continue
           elif test "$i" = mcoll_perf ; then
              for j in `seq 0 9` ; do
                  ext="2.4.$j.nc"
                  bb_ext="bb.2.4.$j.nc"
                  # echo "${LINENO}:--- validating file ${OUT_FILE}.$bb_ext"
                  ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.$bb_ext
                  # echo "${LINENO}:--- ncmpidiff -q ${OUT_FILE}.$ext ${OUT_FILE}.$bb_ext"
                  ${MPIRUN} ${NCMPIDIFF} -q ${OUT_FILE}.$ext ${OUT_FILE}.$bb_ext
              done
              continue
           else
              # echo "--- validating file ${OUT_FILE}.bb.nc"
              ${TESTSEQRUN} ${VALIDATOR} -q ${OUT_FILE}.bb.nc
           fi

           # compare file header only for large file tests
           DIFF_OPT="-q"
           if test "$i" = last_large_var ||
              test "$i" = dim_cdf12 ||
              test "$i" = tst_cdl_hdr_parser ||
              test "$i" = bigrecords ||
              test "$i" = high_dim_var ||
              test "$i" = large_attr ||
              test "$i" = large_coalesce ||
              test "$i" = large_dims_vars_attrs ||
              test "$i" = large_files ||
              test "$i" = large_header ||
              test "$i" = large_reqs ||
              test "$i" = large_var ||
              test "$i" = tst_cdf5_begin ||
              test "$i" = tst_flarge ||
              test "$i" = tst_hash_large_ndims ||
              test "$i" = tst_hash_large_ngattrs ||
              test "$i" = tst_hash_large_nvars ; then
              DIFF_OPT+=" -h"
           fi
           # echo "${LINENO}: --- ncmpidiff $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} $DIFF_OPT $OUT_FILE.nc $OUT_FILE.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${OUT_FILE}.nc4 4
           # Validator does not support nc4
        fi
    done # intra_aggr
    done # mpiio_mode

    if [[ "$i" == *"vard"* ]] ; then
       continue
    fi

    DIFF_OPT="-q"
    if test "$i" = last_large_var ||
       test "$i" = dim_cdf12 ||
       test "$i" = tst_cdl_hdr_parser ||
       test "$i" = bigrecords ||
       test "$i" = high_dim_var ||
       test "$i" = large_attr ||
       test "$i" = large_coalesce ||
       test "$i" = large_dims_vars_attrs ||
       test "$i" = large_files ||
       test "$i" = large_header ||
       test "$i" = large_reqs ||
       test "$i" = large_var ||
       test "$i" = tst_cdf5_begin ||
       test "$i" = tst_flarge ||
       test "$i" = tst_hash_large_ndims ||
       test "$i" = tst_hash_large_ngattrs ||
       test "$i" = tst_hash_large_nvars ; then
       DIFF_OPT+=" -h"
    fi
    if test "$i" = pres_temp_4D_rd || test "$i" = test_inq_format ; then
       continue
    fi
    if test "$i" = put_all_kinds || test "$i" = iput_all_kinds ; then
       for j in 1 2 5; do
          # echo "${LINENO}: --- ncmpidiff $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.mpio.ina.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.mpio.ina.nc$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc$j $OUT_PREFIX.pncio.nc$j
          # echo "--- ncmpidiff $OUT_PREFIX.pncio.nc$j $OUT_PREFIX.pncio.ina.nc$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc$j $OUT_PREFIX.pncio.ina.nc$j
       done
    elif test "$i" = tst_pthread ; then
       for j in `seq 0 ${NTHREADS}` ; do
          # echo "${LINENO}: --- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.mpio.ina.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.nc.$j $OUT_PREFIX.pncio.nc.$j
          # echo "--- ncmpidiff $OUT_PREFIX.pncio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.nc.$j $OUT_PREFIX.pncio.ina.nc.$j
       done
    elif test "$i" = mcoll_perf ; then
       for j in `seq 0 9` ; do
          ext="2.4.$j.nc"
          # echo "${LINENO}: --- ncmpidiff $OUT_PREFIX.mpio.$ext $OUT_PREFIX.mpio.ina.$ext ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.$ext $OUT_PREFIX.mpio.ina.$ext
          # echo "--- ncmpidiff $OUT_PREFIX.mpio.$ext $OUT_PREFIX.pncio.$ext ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.mpio.$ext $OUT_PREFIX.pncio.$ext
          # echo "--- ncmpidiff $OUT_PREFIX.pncio.$ext $OUT_PREFIX.pncio.ina.$ext ---"
          $MPIRUN $NCMPIDIFF $DIFF_OPT $OUT_PREFIX.pncio.$ext $OUT_PREFIX.pncio.ina.$ext
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

    if test "x$i" = xpres_temp_4D_wr ; then
       rm -f ${OUTDIR}/pres_temp_4D*.nc*
    else
       rm -f ${OUTDIR}/$i*nc*
    fi
done # check_PROGRAMS

