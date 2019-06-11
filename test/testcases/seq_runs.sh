#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

${TESTSEQRUN} ./tst_version

${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc1
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc2
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc5

${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds.nc
${TESTSEQRUN} ${VALIDATOR} -q  ${TESTOUTDIR}/iput_all_kinds.nc1
${TESTSEQRUN} ${VALIDATOR} -q  ${TESTOUTDIR}/iput_all_kinds.nc2
${TESTSEQRUN} ${VALIDATOR} -q  ${TESTOUTDIR}/iput_all_kinds.nc5

NCMPIGEN=../../src/utils/ncmpigen/ncmpigen
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# remove the file system type prefix name if there is any.
OUT_PATH=`echo "$TESTOUTDIR" | cut -d: -f2-`

rm -f ${OUT_PATH}/testfile.nc ${OUT_PATH}/redef1.nc
${TESTSEQRUN} ${NCMPIGEN} -v 5 -o ${TESTOUTDIR}/redef1.nc ${srcdir}/redef-good.ncdump
${TESTSEQRUN} ./redef1 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
# diff -q ${OUT_PATH}/testfile.nc ${OUT_PATH}/redef1.nc

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

# echo ""

if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
   # echo "---- testing burst buffering"
   # Run using burst buffer driver
   export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
   ${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds.bb.nc
   ${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds.bb.nc
   unset PNETCDF_HINTS

   # Compare
   for i in 1 2 5 ; do
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/put_all_kinds.nc$i ${TESTOUTDIR}/put_all_kinds.bb.nc$i
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/iput_all_kinds.nc$i ${TESTOUTDIR}/iput_all_kinds.bb.nc$i
   done
fi

# echo ""

if test "${ENABLE_THREAD_SAFE}" = 1 ; then
   # echo "---- testing thread safety"
   for j in 0 1 ; do
       export PNETCDF_SAFE_MODE=$j
       # echo "---- set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"

       ${TESTSEQRUN} ./tst_pthread ${TESTOUTDIR}/tst_pthread.nc
       for i in 0 1 2 3 4 5 ; do
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_pthread.nc.$i
       done
   done
fi
