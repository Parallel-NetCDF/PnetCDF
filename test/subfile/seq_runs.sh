#!/bin/sh

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in $TESTPROGRAMS; do { \
        ${TESTSEQRUN} ./$i         -f ${TESTOUTDIR}/testfile.nc -s 2 ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc ; \
} ; done ; } ; done

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc.subfile_0.nc
