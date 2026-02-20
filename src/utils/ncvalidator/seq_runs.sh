#!/bin/sh
#
# Copyright (C) 2017, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${ENULLPAD_FILES} ; do
    ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_ENULLPAD
done

for i in ${EMAXVARS_FILES} ; do
    ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EMAXVARS
done

for i in ${EUNLIMIT_FILES} ; do
    ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EUNLIMIT
done

for i in ${ENOTNC_FILES} ; do
    ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_ENOTNC
done

for i in ${EVARSIZE_FILES} ; do
    ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EVARSIZE
done

