#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

# Run without DataWarp driver
${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds.nc
${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds.nc

# Run using DataWarp driver
export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable"
${TESTSEQRUN} ./put_all_kinds ${TESTOUTDIR}/put_all_kinds_dw.nc
${TESTSEQRUN} ./iput_all_kinds ${TESTOUTDIR}/iput_all_kinds_dw.nc
unset PNETCDF_HINTS

# Compare
for i in 1 2 5 ; do
# echo "${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf$i ${TESTOUTDIR}/put_all_kinds_dw.nc.cdf$i"
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/put_all_kinds.nc.cdf$i ${TESTOUTDIR}/put_all_kinds_dw.nc.cdf$i
# diff -q ${TESTOUTDIR}/put_all_kinds.nc.cdf$i ${TESTOUTDIR}/put_all_kinds_dw.nc.cdf$i
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf$i ${TESTOUTDIR}/iput_all_kinds_dw.nc.cdf$i
# diff -q ${TESTOUTDIR}/iput_all_kinds.nc.cdf$i ${TESTOUTDIR}/iput_all_kinds_dw.nc.cdf$i
done


