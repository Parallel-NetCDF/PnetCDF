#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

${TESTSEQRUN} ./tst_io ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_io1.nc

if [ -n "${TESTBB}" ]; then
    export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable"
    ${TESTSEQRUN} ./tst_io ${TESTOUTDIR}
    unset PNETCDF_HINTS
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_io1.nc
fi
