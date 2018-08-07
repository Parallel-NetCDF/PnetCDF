#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf1
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf2
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf5

${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf1
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf2
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf5

NCMPIGEN=../../src/utils/ncmpigen/ncmpigen
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# remove the file system type prefix name if there is any.
OUT_PATH=`echo "$TESTOUTDIR" | cut -d: -f2-`

rm -f ${OUT_PATH}/testfile.nc ${OUT_PATH}/redef1.nc
${TESTSEQRUN} ${NCMPIGEN} -v 2 -o ${TESTOUTDIR}/redef1.nc ${srcdir}/redef-good.ncdump
${TESTSEQRUN} ./redef1 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
# diff -q ${OUT_PATH}/testfile.nc ${OUT_PATH}/redef1.nc

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

if [ -n "${TESTBB}" ]; then
   # Run using burst buffer driver
   export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
   ${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds_bb.nc
   ${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds_bb.nc
   unset PNETCDF_HINTS

   # Compare
   for i in 1 2 5 ; do
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf$i ${TESTOUTDIR}/put_all_kinds_bb.nc.cdf$i
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf$i ${TESTOUTDIR}/iput_all_kinds_bb.nc.cdf$i
   done
fi

if [ -n "${TEST_THREAD_SAFE}" ]; then
   for j in 0 1 ; do
       export PNETCDF_SAFE_MODE=$j
       ${TESTSEQRUN} ./tst_pthread ${TESTOUTDIR}/tst_pthread.nc
       for i in 0 1 2 3 4 5 ; do
           ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_pthread.nc.$i
       done
   done
fi
