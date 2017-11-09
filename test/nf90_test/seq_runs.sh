#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

rm -f ${TESTOUTDIR}/scratch.nc ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nf90_test -c    -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test       -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}      -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nf90_test -c -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -2    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}      -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nf90_test -c -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ./nf90_test -5    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}      -q ${TESTOUTDIR}/test.nc


