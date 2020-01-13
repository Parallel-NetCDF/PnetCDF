#!/bin/sh
#
# Copyright (C) 2017, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

for i in ${ENULLPAD_FILES} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_ENULLPAD
    done
done

for i in ${EMAXVARS_FILES} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EMAXVARS
    done
done

for i in ${EUNLIMIT_FILES} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EUNLIMIT
    done
done

for i in ${ENOTNC_FILES} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_ENOTNC
    done
done

for i in ${EVARSIZE_FILES} ; do
    for j in 0 1 ; do
        export PNETCDF_SAFE_MODE=$j
        # echo "set PNETCDF_SAFE_MODE ${PNETCDF_SAFE_MODE}"
        ${TESTSEQRUN} ./tst_open ${srcdir}/$i NC_EVARSIZE
    done
done

