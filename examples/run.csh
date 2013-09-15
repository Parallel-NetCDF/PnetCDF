#!/bin/csh
#
# Copyright (C) 2013, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@


set MPI_RUN="mpiexec -n 4"
set FS = /orangefs/wkliao/testfile.nc

set EXE_DIR=${PWD}

# run example programs in directory C
foreach exe (collective_write \
             nonblocking_write)
    echo "---- ${EXE_DIR}/C/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/C/$exe 10 $FS
    echo""
end 

foreach exe (column_wise \
             block_cyclic \
             flexible_api \
             get_info \
             hints \
             mput \
             put_vara \
             put_varn_float \
             put_varn_int)
    echo "---- ${EXE_DIR}/C/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/C/$exe $FS
    echo""
end 


# run example programs in directory F77
foreach exe (nonblocking_write)
    echo "---- ${EXE_DIR}/F77/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/F77/$exe 10 $FS
    echo""
end 

foreach exe (block_cyclic \
             column_wise \
             flexible_api \
             get_info \
             hints \
             put_vara \
             put_varn_int \
             put_varn_real)
    echo "---- ${EXE_DIR}/F77/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/F77/$exe $FS
    echo""
end 


# run example programs in directory F90
foreach exe (nonblocking_write)
    echo "---- ${EXE_DIR}/F90/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/F90/$exe 10 $FS
    echo""
end 

foreach exe (block_cyclic \
             column_wise \
             flexible_api \
             get_info \
             hints \
             put_vara \
             put_varn_int \
             put_varn_real)
    echo "---- ${EXE_DIR}/F90/$exe ----"
    ${MPI_RUN} ${EXE_DIR}/F90/$exe $FS
    echo""
end 


