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

# remove file system type prefix if there is any
OUTDIR=`echo "$TESTOUTDIR" | cut -d: -f2-`

# let NTHREADS=$1*6-1
NTHREADS=`expr $1 \* 6 - 1`

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test ${PNETCDF_DEBUG} = 1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    for j in ${safe_modes} ; do
    for intra_aggr in 0 1 ; do
        if test "$j" = 1 ; then # test only in safe mode
           export PNETCDF_HINTS="romio_no_indep_rw=true"
        else
           export PNETCDF_HINTS=
        fi
        if test "$intra_aggr" = 1 ; then
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_num_aggrs_per_node=2"
        fi
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        if test "$i" = tst_version ; then
           ${MPIRUN} ./tst_version
           continue
        fi

        if test "$i" = tst_pthread ; then
           # each MPI process created 6 threads
           ${MPIRUN} ./tst_pthread ${TESTOUTDIR}/tst_pthread.nc
           for k in `seq 0 ${NTHREADS}` ; do
               ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_pthread.nc.$k
               rm -f ${OUTDIR}/tst_pthread.nc.$k
           done
           continue
        fi

        ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc

        # put_all_kinds and iput_all_kinds output 3 files
        if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
           for k in 1 2 5 ; do
               # echo "--- validating file ${TESTOUTDIR}/$i.nc$k"
               ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc$k
           done
        else
           # echo "--- validating file ${TESTOUTDIR}/$i.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        fi
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.bb.nc
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # put_all_kinds and iput_all_kinds output 3 files
           if test "$i" = put_all_kinds -o "$i" = iput_all_kinds ; then
              for k in 1 2 5 ; do
                  # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc$k"
                  ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc$k
                  # echo "--- ncmpidiff $i.nc$k $i.bb.nc$k ---"
                  ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc$k ${TESTOUTDIR}/$i.bb.nc$k
              done
              continue
           else
              # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
              ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc
           fi

           # skip ncmpidiff for large file
           if test "$i" = last_large_var ; then
              continue
           fi

           # echo "--- ncmpidiff $i.nc $i.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc ${TESTOUTDIR}/$i.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc4 4
           # Validator does not support nc4
        fi
    done
    done
    rm -f ${OUTDIR}/$i.nc*
    rm -f ${OUTDIR}/$i.bb.nc*
done

