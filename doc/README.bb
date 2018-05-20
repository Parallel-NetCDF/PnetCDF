#
# Copyright (C) 2017, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$

-----------------------------------------------------------------------------
 Using burst buffer as Burst Buffers in PnetCDF
-----------------------------------------------------------------------------

burst buffer driver is a I/O driver in PnetCDF that implements a log-based I/O
aggregation for write related I/O requests that is designed to work on Cray
burst buffer.

-----------------------------------------------------------------------------
 Build PnetCDF with burst buffer driver
-----------------------------------------------------------------------------

To build PnetCDF with burst buffer driver support, simply set "--enable-burstbuffer"
option at configure time:

./configure --prefix=/path/to/install --enable-burstbuffer

-----------------------------------------------------------------------------
 Running with burst buffer Driver
-----------------------------------------------------------------------------

The burst buffer driver is enable by setting file hints on file creation/opening.
To enable burst buffer driver, set the hint "nc_bb" to enable.

MPI_Info_set(info, "nc_bb", "enable");

The hint can also be set using environment variable PNETCDF_HINTS.

export PNETCDF_HINTS="nc_bb=enable"

-----------------------------------------------------------------------------
 Using PnetCDF with burst buffer driver
-----------------------------------------------------------------------------

The burst buffer can be configured using hints. Here's a list of supported hints:

Hint                    Values          Default  Description
----                    ------          -------  -----------
nc_bb                   enable/disable  disable  Whether burst buffer driver is enabled.
nc_bb_dirname           <Valid POSIX    ./       Directory where log file will be
                        Directory>               stored. This is usually set to the
                                                 path where burst buffer is mounted.
nc_bb_del_on_close      enable/disable  enable   Whether logfile should be deleted
                                                 after closing the NetCDF file. It
                                                 can be disabled when the scheduler
                                                 will clean up the burst buffer
                                                 automatically after the job is
                                                 completed.
nc_bb_flush_buffer_size <integer>       0        Amount of memory that can be used
                                                 to flush the log. The unit is in
                                                 bytes. 0 means unlimited. Any write
                                                 request that is larger than the
                                                 buffer size will not be buffered,
                                                 instead, it will be written to PFS
                                                 directly.

-----------------------------------------------------------------------------
 Submitting Job that Enables burst buffer Driver
-----------------------------------------------------------------------------

We show an example script for enabling burst buffer driver on Cori at NERSC

#!/bin/bash 
#SBATCH -p debug 
#SBATCH -N 1 
#SBATCH -C haswell 
#SBATCH -t 00:10:00 
#SBATCH -o output.txt 
#DW jobdw capacity=1289GiB access_mode=private type=scratch pool=sm_pool 
export PNETCDF_HINTS="nc_bb=enable;nc_bb_del_on_close=disable;nc_bb_dirname=${BB_JOB_PRIVATE}" 
srun -n 1 ./myapplication 

-----------------------------------------------------------------------------
 Known Problems
-----------------------------------------------------------------------------
While we design the burst buffer driver to be as transparent as possible. There are
some behaviors that can change when the burst buffer driver is used. Here's a list
of different behaviors:

1. Log buffering delays actual file write to replay time. If there are errors
caused by put operations, it will be hide by logging until the log is replayed.

   