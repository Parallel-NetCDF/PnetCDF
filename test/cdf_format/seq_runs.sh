#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

${TESTSEQRUN} ./test_inq_format ${srcdir}

# the followings check files with corrupted header
${TESTSEQRUN} ./tst_open_cdf5 ${srcdir}/bad_begin.nc5
${TESTSEQRUN} ./tst_corrupt ${srcdir}

# echo ""

if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
    # echo "---- testing burst buffering"

    export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"

    ${TESTSEQRUN} ./test_inq_format ${srcdir}

    # the followings check files with corrupted header
    ${TESTSEQRUN} ./tst_open_cdf5 ${srcdir}/bad_begin.nc5
    ${TESTSEQRUN} ./tst_corrupt ${srcdir}
fi
