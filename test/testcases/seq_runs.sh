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
diff -q ${OUT_PATH}/testfile.nc ${OUT_PATH}/redef1.nc

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

if [ -n "${TESTDW}" ]; then
   # Run using DataWarp driver
   export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable"
   ${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds_dw.nc
   ${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds_dw.nc
   unset PNETCDF_HINTS

   # Compare
   for i in 1 2 5 ; do
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf$i ${TESTOUTDIR}/put_all_kinds_dw.nc.cdf$i
       ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf$i ${TESTOUTDIR}/iput_all_kinds_dw.nc.cdf$i
   done
fi

