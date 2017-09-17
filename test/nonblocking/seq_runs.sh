#!/bin/sh

# "set -x" expands variables and prints a little + sign before the line
set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in ${TESTPROGRAMS}; do ( \
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/testfile.nc ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc ; \
) ; done ; } ; done


${TESTSEQRUN} ./mcoll_perf ${TESTOUTDIR}/testfile
for j in `seq 0 9` ; do { \
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.2.4.$j.nc ; \
} ; done
