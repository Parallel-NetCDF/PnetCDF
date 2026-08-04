#!/bin/bash  -l

# set a file with progressive striping (5 components)
# 1st component is from offset   0 to 4MB       using striping count 4
# 2nd component is from offset   4MB to   1GB   using striping count 8
# 2nd component is from offset   1GB to  10GB   using striping count 16
# 3rd component is from offset  10GB to 100GB   using striping count 32
# 3rd component is from offset 100GB to ~       using striping count -1, i.e. all
# lfs setstripe -E 4m -c 4 -E 1g -c 8 -E 10g -c 16 -E 100g -c 32 -E eof -c -1  /eagle/radix-io/wkliao/pfl
# cp /eagle/radix-io/wkliao/map_i_case_fill_1344p.nc /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc

# set -x

INSTALL_DIR=$HOME/PnetCDF/Github

OPTS=-O2

LUSTRE_DIR=/eagle/radix-io/wkliao

# export LD_LIBRARY_PATH="$HOME/PnetCDF/Github/lib:${LD_LIBRARY_PATH}"

EXE_FILE=tst_pfl

SRC_DIR=.

mpicc -Wall $OPTS -o ${EXE_FILE} $SRC_DIR/$EXE_FILE.c \
                  -I${INSTALL_DIR}/include \
                  -L${INSTALL_DIR}/lib -lpnetcdf -llustreapi

exit

./${EXE_FILE} /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32
