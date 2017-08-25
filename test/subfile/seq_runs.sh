#!/bin/sh

for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$i -f $TESTOUTDIR/testfile.nc -s 2 \
        ; ) ; done
