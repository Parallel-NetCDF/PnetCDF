#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
# set -e
# Cannot use "set -e" here, as the tests here all return errors.

VALIDATOR=ncmpivalid
if [ -z "${VALIDATOR}" ]; then
   echo "Failed: variable VALIDATOR id not defined"
   exit 1
fi
if [ ! -f ${VALIDATOR} ]; then
   echo "Failed: file ${VALIDATOR} does not exit"
   exit 1
fi

for j in 0 1 ; do
    export PNETCDF_SAFE_MODE=$j
    for i in ${BAD_FILES} ; do
        if [ ! -f $i ]; then
           echo "Failed: input test file $i does not exit"
           exit 1
        fi
        ${TESTSEQRUN} ./${VALIDATOR} -q $i
        # capture exit status of VALIDATOR command
        if [ $? -eq 0 ]; then
           echo "Failed: ${VALIDATOR} -q $i"
           exit 1
        fi
    done
done
echo "SUCCESS: ${VALIDATOR} to detect files failing to conform CDF formats"

