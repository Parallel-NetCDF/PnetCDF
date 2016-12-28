#!/bin/sh

for i in $TESTPROGRAMS; do ( \
        $TESTSEQRUN ./$$i $TESTOUTDIR/testfile.nc \
        ; ) ; done

