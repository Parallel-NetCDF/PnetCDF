#!/bin/bash  -l

set -e

# https://docs.nersc.gov/jobs/best-practices/#hugepages
# module load craype-hugepages2M

MY_GROUP=`groups | cut -d' ' -f2`
# echo "MY_GROUP=$MY_GROUP"

dry_run=0

PNETCDF_1141_DIR=/global/common/software/$MY_GROUP/$USER/PnetCDF/1.14.1
SRC_DIR=../PnetCDF

CC=cc
CXX=CC
FC=ftn
if test "x$OPENMPI_ROOT" != x ; then
   # 'module load openmpi' will set OPENMPI_ROOT
   CC=mpicc
   CXX=mpicxx
   FC=mpifort
   echo "---- Using OpenMPI to compile ----"
fi
if test "x$LMOD_FAMILY_PRGENV" = "xPrgEnv-gnu" ; then
   FFLAGS="-fallow-argument-mismatch"
   FCFLAGS="$FFLAGS -ffree-form"
fi

# OPTS="-O0 -g -Wall"
OPTS=-O2
LIBS="-lpnetcdf"

if test "x$PNETCDF_DIR" = x ; then
   PNETCDF_DIR=/global/common/software/$MY_GROUP/$USER/PnetCDF/Github
fi
# When module load PrgEnv-intel
# PNETCDF_DIR= /opt/cray/pe/parallel-netcdf/1.12.3.13/intel/2023.2
# When module load PrgEnv-nvidia
# PNETCDF_DIR=/opt/cray/pe/parallel-netcdf/1.12.3.13/nvidia/23.3

echo "LMOD_FAMILY_PRGENV $LMOD_FAMILY_PRGENV"
echo "PNETCDF_DIR=$PNETCDF_DIR"
echo "PNETCDF_1141_DIR=$PNETCDF_1141_DIR"

EXE_EXTS="master 1.14.1"

EXE_FILES="wrf_io block_cyclic check_constant nonblocking_write put_varn_real flash_benchmark_io"

if test "x$1" = xclean ; then
   for exe in $EXE_FILES ; do
   for ext in $EXE_EXTS ; do
      echo "rm -f $exe.$ext"
      if test "x$dry_run" = x0 ; then
         rm -f $exe.$ext
      fi
   done
   done
   exit
fi

for exe in $EXE_FILES ; do
for ext in $EXE_EXTS ; do

   if test "x$ext" = xmaster ; then
      install_dir=$PNETCDF_DIR
   else
      install_dir=$PNETCDF_1141_DIR
   fi

   echo ""
   echo "==========================================================="
   echo "---- build $exe with PnetCDF $ext"
   echo ""

   has_fortran=`$install_dir/bin/pnetcdf-config --has-fortran`

   INC_DIR="-I $install_dir/include"

   # add -lmpifort when using Intel compilers
   # LIBS+="-lmpifort"

   case "$exe" in
      wrf_io)
          SRC=$SRC_DIR/benchmarks/WRF-IO/${exe}.c
          MPI_COMPILER=$CC
          COMPILE_FLAGS="$OPTS $CFLAGS $INC_DIR"
          ;;
      block_cyclic)
          SRC=$SRC_DIR/examples/CXX/${exe}.cpp
          MPI_COMPILER=$CXX
          COMPILE_FLAGS="$OPTS $CXXFLAGS $INC_DIR"
          ;;
      check_constant)
          SRC=check_constant.f90
          MPI_COMPILER=$FC
          COMPILE_FLAGS="$OPTS $FCFLAGS $INC_DIR"
          ;;
      nonblocking_write)
          SRC="$SRC_DIR/examples/F77/utils.F90 $SRC_DIR/examples/F77/${exe}.f"
          MPI_COMPILER=$FC
          COMPILE_FLAGS="$OPTS $FFLAGS $INC_DIR"
          ;;
      put_varn_real)
          SRC="$SRC_DIR/examples/F90/utils.F90 $SRC_DIR/examples/F90/${exe}.f90"
          MPI_COMPILER=$FC
          COMPILE_FLAGS="$OPTS $FCFLAGS $INC_DIR"
          ;;
      flash_benchmark_io)
          SRC="$SRC_DIR/benchmarks/FLASH-IO/get_mfluid_property.F90 \
               $SRC_DIR/benchmarks/FLASH-IO/flash_release.F90 \
               $SRC_DIR/benchmarks/FLASH-IO/flash_benchmark_io.F90 \
               $SRC_DIR/benchmarks/FLASH-IO/checkpoint_ncmpi_parallel.F90 \
               $SRC_DIR/benchmarks/FLASH-IO/plotfile_ncmpi_parallel.F90"
          MPI_COMPILER=$FC
          COMPILE_FLAGS="$OPTS $FCFLAGS $INC_DIR -I$SRC_DIR/benchmarks/FLASH-IO"
          COMPILE_FLAGS+=" -DN_DIM=3 -DMAXBLOCKS=100 -DIONMAX=13"
          ;;
      *)
          ;;
   esac

   CMD="$MPI_COMPILER $COMPILE_FLAGS $SRC -o ${exe}.$ext -L$install_dir/lib $LIBS"
   echo "CMD=$CMD"
   if test "x$dry_run" = x0 ; then
      $CMD
   fi

done # exe
done # ext


