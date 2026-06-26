#!/bin/bash  -l

#SBATCH --constraint=cpu

#SBATCH --qos=regular
#       #SBATCH --qos=debug

#SBATCH -t 00:08:00

#   	#SBATCH --nodes=1
#       #SBATCH --nodes=2
#       #SBATCH --nodes=4
#       #SBATCH --nodes=6
#       #SBATCH --nodes=8
#       #SBATCH --nodes=16
#       #SBATCH --nodes=32
#SBATCH --nodes=64
#       #SBATCH --nodes=128

#       #SBATCH --job-name=wrf_io.64p
#       #SBATCH --job-name=wrf_io.128
#       #SBATCH --job-name=wrf_io.768p
#       #SBATCH --job-name=wrf_io.1024p
#       #SBATCH --job-name=wrf_io.2048p
#       #SBATCH --job-name=wrf_io.4096p
#SBATCH --job-name=wrf_io.8192p
#       #SBATCH --job-name=wrf_io.16384p

# DO NOT set --ntasks-per-node, as Slurm will assign MPI ranks evenly among
# compute nodes !!!
# When setting --ntasks-per-node=64, and running NP=224 on 4 nodes will get
# warning message: srun: warning: can't honor --ntasks-per-node set to 64 which
# doesn't match the requested tasks 224 with the number of requested nodes 4.
# Ignoring --ntasks-per-node. This also gets you 56 MPI processes per node.
# Another example. Seeting NP=226 will get you 57 + 57 + 56 + 56 on 4 nodes.
#       #SBATCH --ntasks-per-node=64

#SBATCH -o qout.%x.%j
#SBATCH -e qout.%x.%j

#       #SBATCH --mail-user=
#       #SBATCH --mail-type=ALL
#       #SBATCH --mail-type=BEGIN
#SBATCH --dependency=afterany:54666679
# afternay: This job can begin execution after the specified jobs have terminated.
#------------------------------------------------------------------------#
ulimit -c unlimited
# module unload darshan

cd $PWD

MY_GROUP=`groups | cut -d' ' -f2`
# echo "MY_GROUP=$MY_GROUP"

# https://docs.nersc.gov/jobs/best-practices/#hugepages
# module load craype-hugepages2M

DRY_RUN=0
if test "x$SLURM_JOB_QOS" = x ; then
   DRY_RUN=1
fi
if test $DRY_RUN = 1 ; then
   SLURM_JOB_NAME=nonblocking_write
   SLURM_JOB_NAME=mpio_aggre
   SLURM_JOB_NAME=flash_io
   SLURM_JOB_NAME=aggre
   SLURM_JOB_NAME=wrf_io

   SLURM_JOB_NUM_NODES=64
   SLURM_JOB_NUM_NODES=32
   SLURM_JOB_NUM_NODES=4
   SLURM_JOB_NUM_NODES=6
   SLURM_JOB_NUM_NODES=1
   SLURM_JOB_NUM_NODES=16
   SLURM_JOB_NUM_NODES=8
fi

if test "x$SLURM_NTASKS_PER_NODE" = x ; then
   SLURM_NTASKS_PER_NODE=64
   SLURM_NTASKS_PER_NODE=128
fi
NP=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE))

# !!! When running -y 1300 -x 1900, using 2 nodes, 64 processes per node, if
# FI_MR_CACHE_MONITOR is set to kdreg2 and FI_CXI_RX_MATCH_MODE to software
# will cause program hang. Later iterations of ntimes in the two-phase I/O,
# some processes do not complete.

# Sometimes kdreg2 is necessary for F case running 169 nodes, 128 processes per node
# export FI_MR_CACHE_MONITOR=kdreg2
# export FI_UNIVERSE_SIZE=$NP
# export FI_CXI_DEFAULT_CQ_SIZE=524288
# export FI_CXI_RX_MATCH_MODE=software
# setting FI_CXI_RX_MATCH_MODE to software is critical on Perlmutter
# See https://github.com/E3SM-Project/E3SM/issues/5057
# FI_UNIVERSE_SIZE: No. ranks an endpoint would communicate with (default: 256).
# export MPICH_COLL_SYNC=MPI_Bcast
# export MPICH_OFI_NIC_POLICY=NUMA

echo "------------------------------------------------------"
echo "---- Running on Perlmutter CPU nodes ----"
echo "---- SLURM_CLUSTER_NAME      = $SLURM_CLUSTER_NAME"
echo "---- SLURM_JOB_QOS           = $SLURM_JOB_QOS"
echo "---- SLURM_JOB_PARTITION     = $SLURM_JOB_PARTITION"
echo "---- SLURM_JOB_NAME          = $SLURM_JOB_NAME"
echo "---- SBATCH_CONSTRAINT       = $SBATCH_CONSTRAINT"
echo "---- SLURM_JOB_NODELIST      = $SLURM_JOB_NODELIST"
echo "---- SLURM_JOB_NUM_NODES     = $SLURM_JOB_NUM_NODES"
echo "---- SLURM_NTASKS_PER_NODE   = $SLURM_NTASKS_PER_NODE"
echo "---- SLURM_JOB_ID            = $SLURM_JOB_ID"
echo "---- SLURM out/err file      = qout.$SLURM_JOB_NAME.$SLURM_JOB_ID"
echo ""
echo "ENV explicitly set:"
echo "---- FI_MR_CACHE_MONITOR     = $FI_MR_CACHE_MONITOR"
echo "---- FI_UNIVERSE_SIZE        = $FI_UNIVERSE_SIZE"
echo "---- FI_CXI_DEFAULT_CQ_SIZE  = $FI_CXI_DEFAULT_CQ_SIZE"
echo "---- FI_CXI_RX_MATCH_MODE    = $FI_CXI_RX_MATCH_MODE"
echo "---- MPICH_COLL_SYNC         = $MPICH_COLL_SYNC"
echo "---- MPICH_OFI_NIC_POLICY    = $MPICH_OFI_NIC_POLICY"
echo "------------------------------------------------------"
echo ""

APP=wrf_io

# For fast executable loading on Cori and Perlmutter
# sbcast -v ${EXE_DIR}/${EXE_FILE} ${EXE}
for ext in master 1.14.1 ; do
   EXE_FILE=${APP}.$ext
   if [ ! -f $EXE_FILE ]; then
      echo "Error: $EXE_FILE executable file does not exist."
      exit 1
   fi
   EXE=/tmp/${USER}_${EXE_FILE}
   SBCAST_CMD="sbcast ${EXE_FILE} ${EXE}"
   echo "SBCAST_CMD=$SBCAST_CMD"
   if test $DRY_RUN = 0 ; then
      $SBCAST_CMD
   fi
done
echo ""

for ext in Github 1.14.1 ; do
   lib_dir=/global/common/software/$MY_GROUP/$USER/PnetCDF/$ext
   if [ ! -d $lib_dir ]; then
      echo "Error: $lib_dir does not exist."
      exit 1
   fi
done

saved_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

striping_unit=1048576
striping_factor=$SLURM_JOB_NUM_NODES

INTRA_NODES="0 32"
INTRA_NODES="0"
NUM_RECS=1

if test $SLURM_JOB_NUM_NODES = 8 ; then
   longitude=2600
   latitude=3800
elif test $SLURM_JOB_NUM_NODES = 16 ; then
   longitude=2600
   latitude=7600
elif test $SLURM_JOB_NUM_NODES = 32 ; then
   longitude=5200
   latitude=7600
elif test $SLURM_JOB_NUM_NODES = 64 ; then
   longitude=5200
   latitude=15200
elif test $SLURM_JOB_NUM_NODES = 128 ; then
   longitude=10400
   latitude=15200
fi

IO_DIR=${SCRATCH}/WRF-IO
CDL_FILE="${IO_DIR}/wrf_header.txt"
OUT_FILE="${IO_DIR}/${APP}.${longitude}x${latitude}.nc"

DIR_STRIPING=`lfs getstripe -d ${IO_DIR}`
# stripe_count:  1 stripe_size:   1048576 pattern:       raid0 stripe_offset: -1 pool:          original
# stripe_count:  8 stripe_size:   1048576 pattern:       raid0,overstriped stripe_offset: -1
lfs_stripe_count=`echo ${DIR_STRIPING} | cut -d' ' -f2`
lfs_stripe_size=`echo ${DIR_STRIPING} | cut -d' ' -f4`
lfs_pattern=`echo ${DIR_STRIPING} | cut -d' ' -f6`
lfs_stripe_offset=`echo ${DIR_STRIPING} | cut -d' ' -f8`
lfs_pool=`echo ${DIR_STRIPING} | cut -d' ' -f10`

IN_FILE_SUFFIX="${longitude}x${latitude}.nc"
FILE_ID=0

DO_READ=true
DO_READ=false
DO_WRITE=false
DO_WRITE=true

NTIMES=3
for ntime in $(seq 1 ${NTIMES}) ; do

for intranode in ${INTRA_NODES} ; do

for ext in master 1.14.1 ; do

   EXE=/tmp/${USER}_${APP}.$ext

   striping_factor=$SLURM_JOB_NUM_NODES
   cb_nodes=${striping_factor}

   hints="striping_factor=${striping_factor};"
   hints+="striping_unit=${striping_unit};"

   if test "x$ext" = xmaster ; then
      lib_dir=/global/common/software/$MY_GROUP/$USER/PnetCDF/Github
   else
      lib_dir=/global/common/software/$MY_GROUP/$USER/PnetCDF/$ext
      hints+="cb_nodes=${cb_nodes};"
   fi

   if test "x$intranode" != x0 ; then
      hints+="nc_num_aggrs_per_node=${intranode};"
   fi

   export LD_LIBRARY_PATH="$lib_dir/lib:$saved_LD_LIBRARY_PATH"
   echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
   echo ""

   export PNETCDF_HINTS=$hints

IN_FILE_SUFFIX="${longitude}x${latitude}.nc"
IN_FILE="$IO_DIR/$FILE_ID.$IN_FILE_SUFFIX"
((FILE_ID++))

echo ""
echo "-------------------------------------------------------------------------"
echo "    EXE_FILE       = ${EXE_FILE}"
echo "    No. nodes      = $SLURM_JOB_NUM_NODES"
echo "    No. MPI ranks  = $NP"
echo "    exe            = ${EXE}"
echo "    ntime          = $ntime out of ${NTIMES}"
echo "    CDL_FILE       = $CDL_FILE"
if test "x${DO_WRITE}" = xtrue ; then
echo "    OUT_FILE       = $OUT_FILE"
fi
if test "x${DO_READ}" = xtrue ; then
echo "    IN_FILE        = $IN_FILE"
fi
echo ""
echo "    I/O folder stripe pool    = $lfs_pool"
echo "    I/O folder stripe pattern = $lfs_pattern"
echo "    I/O folder stripe count   = $lfs_stripe_count"
echo "    I/O folder stripe size    = $lfs_stripe_size"
echo "    I/O folder stripe offset  = $lfs_stripe_offset"
echo ""
echo "    striping_unit         = $striping_unit"
echo "    striping_factor       = $striping_factor"
echo "    start_iodevice        = $start_iodevice"
echo "    overstriping          = $overstriping"
echo "    ncb_ratio             = $ncb_ratio"
echo "    cb_nodes              = $cb_nodes"
echo "    cb_buffer_size        = $bsize"
echo "    nc_num_aggrs_per_node = $intranode"
echo "    PNETCDF_HINTS         = $PNETCDF_HINTS"
echo "-------------------------------------------------------------------------"
echo ""

date

CMD_OPTS="-n $NUM_RECS"

if test "x$DO_WRITE" = xtrue ; then
   OUT_OPTS="-y $longitude -x $latitude -i $CDL_FILE -w $OUT_FILE"
else
   OUT_OPTS=
fi

if test "x$DO_READ" = xtrue ; then
   if test "x$DO_WRITE" = xtrue ; then
      IN_OPTS="-r $OUT_FILE"
   else
      IN_OPTS="-r $IN_FILE"
   fi
else
   IN_OPTS=
fi

if test "x$OUT_OPTS" != x && test "x$INT_OPTS" = x ; then
   # !!! delete the output file now is faster than clobber it later in ncmpi_create() !!!
   if test $DRY_RUN = 0 ; then
      rm -rf ${OUT_FILE}
   else
      echo "rm -rf ${OUT_FILE}"
   fi
fi
echo ""

CMD="srun -n $NP ${EXE} $CMD_OPTS $OUT_OPTS $IN_OPTS"
echo "CMD=$CMD"
if test $DRY_RUN = 0 ; then
   $CMD
fi
echo ""

echo "---------------------------------------------------------------"
if test $DRY_RUN = 0 ; then
   if test "x$OUT_OPTS" != x ; then
      ls -lR ${OUT_FILE}
      ls -lRh ${OUT_FILE}
      # lfs getstripe -p -L -c -S -i ${OUT_FILE}
      lfs getstripe ${OUT_FILE}
   fi
   if test "x$IN_OPTS" != x ; then
      ls -lR ${IN_FILE}
      ls -lRh ${IN_FILE}
      lfs getstripe -p -L -c -S -i ${IN_FILE}
   fi
else
   if test "x$OUT_OPTS" != x ; then
      echo "ls -lR ${OUT_FILE}"
      echo "ls -lRh ${OUT_FILE}"
   fi
   if test "x$IN_OPTS" != x ; then
      echo "ls -lR ${IN_FILE}"
      echo "ls -lRh ${IN_FILE}"
   fi
fi
echo ""

echo "========================================================================="

done  # loop ext
done  # loop intranode
done  # loop ntime
date
echo ""

