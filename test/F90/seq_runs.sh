#!/bin/sh

set -e

${TESTSEQRUN} ./test_intent  ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_f90      ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./f90tst_vars  ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_types2   ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_f90_cdf5 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./f90tst_vars2 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./f90tst_vars3 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./f90tst_vars4 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_io       ${TESTOUTDIR}

if test "x${RUN_LARGE_FILE_TEST}" = xyes ; then
   ${TESTSEQRUN} ./tst_flarge   ${TESTOUTDIR}/testfile.nc
fi
