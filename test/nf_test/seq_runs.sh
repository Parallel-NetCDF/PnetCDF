#!/bin/sh

set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

rm -f ${TESTOUTDIR}/scratch.nc ${TESTOUTDIR}/test.nc
${TESTSEQRUN} ./nf_test -c    -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf_test       -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nf_test -c -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf_test -2    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nf_test -c -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf_test -5    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc


