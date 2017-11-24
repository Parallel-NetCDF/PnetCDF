#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

# tst_nofill.c creates two files: tst_nofill.nc.fill and tst_nofill.nc.nofill
${TESTSEQRUN} ./tst_nofill    ${TESTOUTDIR}/tst_nofill.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_nofill.nc.fill
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_nofill.nc.nofill


# disable safe mode, as nc_test already runs slow
export PNETCDF_SAFE_MODE=0

rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nc_test    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nc_test -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nc_test -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

