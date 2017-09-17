#!/bin/sh

set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

export PNETCDF_SAFE_MODE=0
${TESTSEQRUN} ./test_inq_format $srcdir
${TESTSEQRUN} ./cdf_type      ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./dim_cdf12     ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

export PNETCDF_SAFE_MODE=1
${TESTSEQRUN} ./test_inq_format $srcdir
${TESTSEQRUN} ./cdf_type      ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./dim_cdf12     ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

