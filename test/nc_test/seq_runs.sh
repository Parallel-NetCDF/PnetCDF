#!/bin/sh

set -e

rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/testfile.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nc_test    -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test       -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -2 -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -2    -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -5 -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -5    -d ${TESTOUTDIR}
${TESTSEQRUN} ./t_nc             ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_misc         ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_norm         ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_small        ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_names        ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_atts3        ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_atts         ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./tst_nofill       ${TESTOUTDIR}/testfile.nc

