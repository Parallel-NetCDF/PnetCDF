#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in ${TESTPROGRAMS}; do ( \
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/$i.nc ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/$i.nc ; \
) ; done ; } ; done

# tst_nofill.c creates two files: tst_nofill.nc.fill and tst_nofill.nc.nofill
${TESTSEQRUN} ./tst_nofill    ${TESTOUTDIR}/tst_nofill.nc
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_nofill.nc.fill
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/tst_nofill.nc.nofill

rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/testfile.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc

${TESTSEQRUN} ./nc_test    -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test       -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nc_test -2 -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -2    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc

${TESTSEQRUN} ./nc_test -5 -c -d ${TESTOUTDIR}
${TESTSEQRUN} ./nc_test -5    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR}    -q ${TESTOUTDIR}/test.nc

