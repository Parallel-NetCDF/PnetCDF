#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

VALIDATOR=../../src/utils/ncvalidator/ncvalidator

for i in ${BAD_FILES} ; do
    ${TESTSEQRUN} ./${VALIDATOR} -q ${srcdir}/$i
done

