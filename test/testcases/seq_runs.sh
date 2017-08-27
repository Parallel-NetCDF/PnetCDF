#!/bin/sh

set -e

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in $TESTPROGRAMS; do { \
        $TESTSEQRUN ./$i $TESTOUTDIR/testfile.nc ; \
} ; done ; } ; done

NCMPIGEN=../../src/utils/ncmpigen/ncmpigen
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

rm -f ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
${TESTSEQRUN} ${NCMPIGEN} -v 2 -o ${TESTOUTDIR}/redef1.nc ${srcdir}/redef-good.ncdump
${TESTSEQRUN} ./redef1 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
diff -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc

