#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
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

# the followings check files with corrupted header

${TESTSEQRUN} ./tst_open_cdf5 ${srcdir}/bad_begin.nc5

${TESTSEQRUN} ./tst_corrupt ${srcdir}

