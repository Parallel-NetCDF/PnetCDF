#!/bin/sh

set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in ${TESTPROGRAMS}; do ( \
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/testfile.nc ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc ; \
) ; done ; } ; done

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

