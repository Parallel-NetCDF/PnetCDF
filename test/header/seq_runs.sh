#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$i $TESTOUTDIR/testfile.nc ; \
) ; done ; } ; done
