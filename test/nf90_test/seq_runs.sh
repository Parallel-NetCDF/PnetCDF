#!/bin/sh

set -e

rm -f ${TESTOUTDIR}/scratch.nc ${TESTOUTDIR}/test.nc
${TESTSEQRUN} ./nf90_test -c    -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test       -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -c -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -2    -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -c -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -5    -d ${TESTOUTDIR}

