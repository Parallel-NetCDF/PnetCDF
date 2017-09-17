#!/bin/sh

set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in ${TESTPROGRAMS}; do ( \
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/testfile.nc ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc ; \
) ; done ; } ; done

${TESTSEQRUN} ./tst_io ${TESTOUTDIR}

