#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

# disable safe mode, as nf_test already runs slow
export PNETCDF_SAFE_MODE=0

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc


