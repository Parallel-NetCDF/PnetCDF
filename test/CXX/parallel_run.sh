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
# echo "TESTPROGRAMS=${TESTPROGRAMS}"

# echo "PNETCDF_DEBUG = ${PNETCDF_DEBUG}"
if test "x${PNETCDF_DEBUG}" = x1 ; then
   safe_modes="0 1"
else
   safe_modes="0"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${TESTPROGRAMS} ; do
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

        ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc

        # echo "--- validating file ${TESTOUTDIR}/$i.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "---- test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.bb.nc
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # echo "--- validating file ${TESTOUTDIR}/$i.bb.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.nc

           # echo "--- ncmpidiff $i.nc $i.bb.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.nc ${TESTOUTDIR}/$i.bb.nc
        fi

        if test "x${ENABLE_NETCDF4}" = x1 ; then
           # echo "test netCDF-4 feature"
           ${MPIRUN} ./$i ${TESTOUTDIR}/$i.nc
           # Validator does not support NetCDF-4 format
        fi
    done
    done
    done
    rm -f ${OUTDIR}/$i.nc
    rm -f ${OUTDIR}/$i.bb.nc
done

