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

if [ -n "${TESTBB}" ]; then
    export PNETCDF_HINTS="nc_bb=enable;nc_bb_dirname=${TESTOUTDIR};nc_bb_overwrite=enable"
    ${TESTSEQRUN} ./test_inq_format ${srcdir}

    # the followings check files with corrupted header

    ${TESTSEQRUN} ./tst_open_cdf5 ${srcdir}/bad_begin.nc5

    ${TESTSEQRUN} ./tst_corrupt ${srcdir}
    unset PNETCDF_HINTS
fi
