#!/bin/sh

for j in 0 1 ; do { \
    export PNETCDF_SAFE_MODE=$$j ; \
    for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$i -f $TESTOUTDIR/testfile.nc -s 2 \
; ) ; done ; } ; done
