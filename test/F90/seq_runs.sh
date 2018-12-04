#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

${TESTSEQRUN} ./tst_io ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_io1.nc
OUTDIR=$(echo $TESTOUTDIR | cut -d: -f2)
mv ${OUTDIR}/tst_io1.nc ${OUTDIR}/tst_io1.nc0

if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
    export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
    ${TESTSEQRUN} ./tst_io ${TESTOUTDIR}
    unset PNETCDF_HINTS
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_io1.nc

    # echo "--- ncmpidiff tst_io1.nc0 tst_io1.nc ---"
    ${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/tst_io1.nc0 ${TESTOUTDIR}/tst_io1.nc

fi
