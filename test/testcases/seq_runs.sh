#!/bin/sh

set -e

VALIDATOR=../../src/utils/ncmpivalid/ncmpivalid

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in ${TESTPROGRAMS}; do { \
        ${TESTSEQRUN} ./$i            ${TESTOUTDIR}/testfile.nc ; \
        ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc ; \
} ; done ; } ; done

NCMPIGEN=../../src/utils/ncmpigen/ncmpigen
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

rm -f ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
${TESTSEQRUN} ${NCMPIGEN} -v 2 -o ${TESTOUTDIR}/redef1.nc ${srcdir}/redef-good.ncdump
${TESTSEQRUN} ./redef1 ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ${NCMPIDIFF} -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc
diff -q ${TESTOUTDIR}/testfile.nc ${TESTOUTDIR}/redef1.nc

${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/testfile.nc

${TESTSEQRUN} ./tst_open_cdf5 ${srcdir}/bad_begin.nc5

${TESTSEQRUN} ./tst_EMAXDIMS ${srcdir}/bad_ndims.nc1
${TESTSEQRUN} ./tst_EMAXDIMS ${srcdir}/bad_ndims.nc2
${TESTSEQRUN} ./tst_EMAXDIMS ${srcdir}/bad_ndims.nc5

${TESTSEQRUN} ./tst_EMAXATTS ${srcdir}/bad_nattrs.nc1
${TESTSEQRUN} ./tst_EMAXATTS ${srcdir}/bad_nattrs.nc2
${TESTSEQRUN} ./tst_EMAXATTS ${srcdir}/bad_nattrs.nc5

${TESTSEQRUN} ./tst_EBADDIM ${srcdir}/bad_dimid.nc1
${TESTSEQRUN} ./tst_EBADDIM ${srcdir}/bad_dimid.nc2
${TESTSEQRUN} ./tst_EBADDIM ${srcdir}/bad_dimid.nc5


