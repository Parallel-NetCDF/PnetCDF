/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include "pncio.h"

/*----< PNCIO_File_close() >--------------------------------------------------*/
int PNCIO_File_close(PNCIO_File *fh)
{
    int err = NC_NOERR;

    if (fh->is_open) {
        err = close(fh->fd_sys);
        if (err != 0)
            err = ncmpii_error_posix2nc("close");
    }

    if (fh->hints->ranklist != NULL)
        NCI_Free(fh->hints->ranklist);
    if (fh->hints != NULL)
        NCI_Free(fh->hints);
    if (fh->info != MPI_INFO_NULL)
        MPI_Info_free(&(fh->info));
    if (fh->io_buf != NULL)
        NCI_Free(fh->io_buf);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    int i, world_rank;
    double timing[NMEASURES*2], max_t[NMEASURES*2], pread_t;
    MPI_Count max_ntimes, counter[NMEASURES*2], max_c[NMEASURES*2];

    /* print two-phase I/O timing breakdown */
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    for (i=0; i<NMEASURES; i++) {
        timing[i]  = fh->write_timing[i];
        counter[i] = fh->write_counter[i];
        timing[i+NMEASURES]  = fh->read_timing[i];
        counter[i+NMEASURES] = fh->read_counter[i];
    }
    MPI_Reduce(timing,  max_t, NMEASURES*2, MPI_DOUBLE, MPI_MAX, 0, fh->comm);
    MPI_Reduce(counter, max_c, NMEASURES*2, MPI_COUNT,  MPI_MAX, 0, fh->comm);

    pread_t = max_t[NMEASURES+2];
    max_ntimes = max_c[0];

    if (world_rank == 0 && max_ntimes > 0) {
        printf("%s: TWO-PHASE write init %5.2f pwrite %5.2f pread %5.2f post %5.2f hsort %5.2f comm %5.2f collw %5.2f\n",
        __func__, max_t[1], max_t[2], pread_t, max_t[4], max_t[5], max_t[3], max_t[0]);
        printf("%s: TWO-PHASE write ntimes %lld check_hole %lld (total_num %lld nrecv %lld) no check %lld (total_num %lld nrecv %lld)\n",
        __func__, max_c[0], max_c[1], max_c[2], max_c[3], max_c[4], max_c[5], max_c[6]);
    }

    max_ntimes = max_c[NMEASURES];

    if (world_rank == 0 && max_ntimes > 0)
        printf("%s: TWO-PHASE read  init %5.2f pread  %5.2f post %5.2f wait %5.2f collr %5.2f ntimes %lld\n",
        __func__, max_t[NMEASURES+1], max_t[NMEASURES+2], max_t[NMEASURES+4], max_t[NMEASURES+3], max_t[NMEASURES+0], max_ntimes);
#endif

    return err;
}
