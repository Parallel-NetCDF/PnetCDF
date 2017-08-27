#!/bin/sh

set -e

for j in 0 1 ; do { \
export PNETCDF_SAFE_MODE=$$j ; \
for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$i $TESTOUTDIR/testfile.nc ; \
) ; done ; } ; done
