#!/bin/bash
#
# Copyright (C) 2019, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

VALIDATOR=../../src/utils/ncvalidator/ncvalidator
NCMPIDIFF=../../src/utils/ncmpidiff/ncmpidiff

MPIRUN=`echo ${TESTMPIRUN} | ${SED} -e "s/NP/$1/g"`
# echo "MPIRUN = ${MPIRUN}"
# echo "check_PROGRAMS=${check_PROGRAMS}"

if test "x$ENABLE_GIO" = x0 ; then
   IO_MODES="mpiio"
else
   IO_MODES="gio mpiio"
fi

# prevent user environment setting of PNETCDF_HINTS to interfere
unset PNETCDF_HINTS

for i in ${check_PROGRAMS} ; do
    for io_mode in $IO_MODES ; do
        if test "x$io_mode" = xmpiio ; then
           USEMPIO_HINTS="nc_driver=mpiio"
        else
           USEMPIO_HINTS="nc_driver=gio"
        fi
    for intra_aggr in 0 1 ; do
        if test "$intra_aggr" = 1 ; then
           INA_HINTS="nc_num_aggrs_per_node=2"
        else
           INA_HINTS="nc_num_aggrs_per_node=0"
        fi

        PNETCDF_HINTS=
        if test "x$USEMPIO_HINTS" != x ; then
           PNETCDF_HINTS="$USEMPIO_HINTS;$PNETCDF_HINTS"
        fi
        if test "x$INA_HINTS" != x ; then
           PNETCDF_HINTS="$INA_HINTS;$PNETCDF_HINTS"
        fi
        if test "x$MIMIC_LUSTRE" != x1 ; then
           PNETCDF_HINTS="cb_nodes=2;$PNETCDF_HINTS"
        fi

        export PNETCDF_HINTS="$PNETCDF_HINTS"
	    # echo "PNETCDF_HINTS=$PNETCDF_HINTS"

        if test "$i" = open ; then
           ${MPIRUN} ./$i ${srcdir}/arrays.bp
           ${MPIRUN} ./$i ${srcdir}/attributes.bp
           ${MPIRUN} ./$i ${srcdir}/arrays_big.bp
           if test ${ADIOS_VER_GE_1132} = 1 ; then
               ${MPIRUN} ./$i ${srcdir}/attributes_big.bp
           fi
       elif test "$i" = att ; then
          ${MPIRUN} ./$i ${srcdir}/attributes.bp
       else
          ${MPIRUN} ./$i ${srcdir}/arrays.bp
       fi
    done
    done
done

