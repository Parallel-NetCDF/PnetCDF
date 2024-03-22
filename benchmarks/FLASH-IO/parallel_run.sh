#!/bin/sh
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
           export PNETCDF_HINTS="romio_no_indep_rw=true"
        fi
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

        ${MPIRUN} ./$i -q ${TESTOUTDIR}/$i.

        if test $? = 0 ; then
           echo "PASS: F90 parallel run on $1 processes --------------- $i"
        fi

        if test "x${BUILD_BENCHMARKS_IN_PNETCDF}" != x1 ; then
           continue
        fi

        # echo "--- validating file ${TESTOUTDIR}/$i.ncmpi_chk_0000.nc"
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.ncmpi_chk_0000.nc
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.ncmpi_plt_cnt_0000.nc
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.ncmpi_plt_crn_0000.nc
        # echo ""

        if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
           # echo "test burst buffering feature"
           saved_PNETCDF_HINTS=${PNETCDF_HINTS}
           export PNETCDF_HINTS="${PNETCDF_HINTS};nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
           ${MPIRUN} ./$i -q ${TESTOUTDIR}/$i.bb.
           if test $? = 0 ; then
              echo "PASS: F90 parallel run on $1 processes --------------- $i"
           fi
           export PNETCDF_HINTS=${saved_PNETCDF_HINTS}

           # echo "--- validating file ${TESTOUTDIR}/$i.bb.ncmpi_chk_0000.nc"
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.ncmpi_chk_0000.nc
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.ncmpi_plt_cnt_0000.nc
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.bb.ncmpi_plt_crn_0000.nc

           # echo "--- ncmpidiff $i.ncmpi_chk_0000.nc $i.bb.ncmpi_chk_0000.nc ---"
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.ncmpi_chk_0000.nc ${TESTOUTDIR}/$i.bb.ncmpi_chk_0000.nc
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.ncmpi_plt_cnt_0000.nc ${TESTOUTDIR}/$i.bb.ncmpi_plt_cnt_0000.nc
           ${MPIRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/$i.ncmpi_plt_crn_0000.nc ${TESTOUTDIR}/$i.bb.ncmpi_plt_crn_0000.nc
        fi
    done
    rm -f ${OUTDIR}/$i.ncmpi_chk_0000.nc
    rm -f ${OUTDIR}/$i.ncmpi_plt_cnt_0000.nc
    rm -f ${OUTDIR}/$i.ncmpi_plt_crn_0000.nc
    rm -f ${OUTDIR}/$i.bb.ncmpi_chk_0000.nc
    rm -f ${OUTDIR}/$i.bb.ncmpi_plt_cnt_0000.nc
    rm -f ${OUTDIR}/$i.bb.ncmpi_plt_crn_0000.nc
done

