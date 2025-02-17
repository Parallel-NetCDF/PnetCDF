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

#include <pnc_debug.h>
#include <common.h>
#include "pncio.h"

/*----< PNCIO_File_set_view() >-----------------------------------------------*/
/* For PnetCDF, this subroutine is an independent call, because PnetCDF only
 * use the followings.
 *   Argument etype is always MPI_BYTE.
 *   Argument datarep is always "native".
 *   Argument info is always MPI_INFO_NULL.
 */
int PNCIO_File_set_view(PNCIO_File   *fd,
                        MPI_Offset    disp,
                        MPI_Datatype  filetype,
                        MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                        MPI_Count    *offsets,
                        MPI_Count    *lengths
#else
                        MPI_Offset   *offsets,
                        int          *lengths
#endif
)
{
    MPI_Aint i;

assert(filetype == MPI_BYTE);
assert(disp == 0);
fd->filetype = filetype;
fd->disp = 0;

    fd->flat_file.count = npairs;
    fd->flat_file.off   = offsets;
    fd->flat_file.len   = lengths;
    fd->flat_file.idx   = 0;
    fd->flat_file.rem   = (npairs > 0) ? lengths[0] : 0;

    /* Size of fileview must be calculated here, as PnetCDF may coalesce the
     * offset-length pairs in order to make offsets sorted in a monotonically
     * non-decreasing order.
     */
    fd->flat_file.size = 0;
    for (i=0; i<npairs; i++) fd->flat_file.size += lengths[i];

    /* is_contig is redundant to (count <= 1), but convenient */
    fd->flat_file.is_contig = (npairs <= 1);

    return NC_NOERR;
}

