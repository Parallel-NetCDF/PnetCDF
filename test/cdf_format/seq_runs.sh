#!/bin/sh

set -e

${TESTSEQRUN} ./test_inq_format $srcdir
${TESTSEQRUN} ./cdf_type  ${TESTOUTDIR}/testfile.nc
${TESTSEQRUN} ./dim_cdf12 ${TESTOUTDIR}/testfile.nc

