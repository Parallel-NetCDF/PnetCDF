#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

# disable safe mode, as nf_test already runs slow
export PNETCDF_SAFE_MODE=0

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test    -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test -2 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

if test "${ENABLE_NETCDF4}" = 1 ; then
    rm -f ${TESTOUTDIR}/test.nc
    rm -f ${TESTOUTDIR}/scratch.nc
    rm -f ${TESTOUTDIR}/tooth-fairy.nc
    ${TESTSEQRUN} ./nf_test -4 -d ${TESTOUTDIR}
    # Validator does not support nc4
    # ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc
fi

rm -f ${TESTOUTDIR}/test.nc
rm -f ${TESTOUTDIR}/scratch.nc
rm -f ${TESTOUTDIR}/tooth-fairy.nc
${TESTSEQRUN} ./nf_test -5 -d ${TESTOUTDIR}
${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc


if test "x${ENABLE_BURST_BUFFER}" = x1 ; then
    # echo "---- testing burst buffering"

    export PNETCDF_HINTS="nc_burst_buf=enable;nc_burst_buf_dirname=${TESTOUTDIR};nc_burst_buf_overwrite=enable"
    rm -f ${TESTOUTDIR}/test.nc
    rm -f ${TESTOUTDIR}/scratch.nc
    rm -f ${TESTOUTDIR}/tooth-fairy.nc
    ${TESTSEQRUN} ./nf_test    -d ${TESTOUTDIR}
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

    rm -f ${TESTOUTDIR}/test.nc
    rm -f ${TESTOUTDIR}/scratch.nc
    rm -f ${TESTOUTDIR}/tooth-fairy.nc
    ${TESTSEQRUN} ./nf_test -2 -d ${TESTOUTDIR}
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc

    rm -f ${TESTOUTDIR}/test.nc
    rm -f ${TESTOUTDIR}/scratch.nc
    rm -f ${TESTOUTDIR}/tooth-fairy.nc
    ${TESTSEQRUN} ./nf_test -5 -d ${TESTOUTDIR}
    ${TESTSEQRUN} ${VALIDATOR} -q ${TESTOUTDIR}/test.nc
fi
