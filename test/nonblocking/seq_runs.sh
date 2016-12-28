#!/bin/sh

set -e

for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$i $TESTOUTDIR/testfile.nc ; \
) ; done
