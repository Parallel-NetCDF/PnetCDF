#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
# set -e
# Cannot use "set -e" here, as the tests here all return errors.

VALIDATOR=ncvalidator
if [ -z "${VALIDATOR}" ]; then
   echo "Failed: variable VALIDATOR id not defined"
   exit 1
fi
if [ ! -f ${VALIDATOR} ]; then
   echo "Failed: file ${VALIDATOR} does not exit"
   exit 1
fi

for i in ${BAD_FILES} ; do
    if [ ! -f ${srcdir}/$i ]; then
       echo "Failed: input test file ${srcdir}/$i does not exit"
       exit 1
    fi
    ${TESTSEQRUN} ./${VALIDATOR} -q ${srcdir}/$i
    ret=$?
    # capture exit status of VALIDATOR command
    if [ ${ret} -ne 1 ]; then
       echo "Failed: ${VALIDATOR} -q ${srcdir}/$i"
       exit 1
    fi
done
echo "SUCCESS: ${VALIDATOR} to detect files fail to conform CDF formats"

