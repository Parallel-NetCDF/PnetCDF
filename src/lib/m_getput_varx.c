/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include "ncmpidtype.h"
#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

static int
ncmpii_mgetput_varm(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[], /* [num] */
                    int                rw_flag,     /* WRITE_REQ or READ_REQ */
                    int                io_method);  /* COLL_IO or INDEP_IO */

/*----< ncmpi_mput_var() >---------------------------------------------------*/
int
ncmpi_mput_var(int           ncid, 
               int           num, 
               int           varids[],    /* [num] */
               void         *bufs[],      /* [num] */
               MPI_Offset    bufcounts[], /* [num] */
               MPI_Datatype  datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_var1() >--------------------------------------------------*/
int
ncmpi_mput_var1(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_vara() >--------------------------------------------------*/
int
ncmpi_mput_vara(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_vars() >--------------------------------------------------*/
int
ncmpi_mput_vars(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                MPI_Offset* const  strides[],   /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_varm() >--------------------------------------------------*/
int
ncmpi_mput_varm(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                MPI_Offset* const  strides[],   /* [num] */
                MPI_Offset* const  imaps[],     /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mget_var() >---------------------------------------------------*/
int
ncmpi_mget_var(int           ncid, 
               int           num, 
               int           varids[],     /* [num] */
               void         *bufs[],       /* [num] */
               MPI_Offset    bufcounts[],  /* [num] */
               MPI_Datatype  datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_var1() >--------------------------------------------------*/
int
ncmpi_mget_var1(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_vara() >--------------------------------------------------*/
int
ncmpi_mget_vara(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_vars() >--------------------------------------------------*/
int
ncmpi_mget_vars(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                MPI_Offset* const  strides[],    /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_varm() >--------------------------------------------------*/
int
ncmpi_mget_varm(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                MPI_Offset* const  strides[],    /* [num] */
                MPI_Offset* const  imaps[],      /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mput_var_all() >-----------------------------------------------*/
int
ncmpi_mput_var_all(int           ncid, 
                   int           num, 
                   int           varids[],    /* [num] */
                   void         *bufs[],      /* [num] */
                   MPI_Offset    bufcounts[], /* [num] */
                   MPI_Datatype  datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

/*----< ncmpi_mput_var1_all() >----------------------------------------------*/
int
ncmpi_mput_var1_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

/*----< ncmpi_mput_vara_all() >----------------------------------------------*/
int
ncmpi_mput_vara_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

/*----< ncmpi_mput_vars_all() >----------------------------------------------*/
int
ncmpi_mput_vars_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

/*----< ncmpi_mput_varm_all() >----------------------------------------------*/
int
ncmpi_mput_varm_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

/*----< ncmpi_mget_var_all() >-----------------------------------------------*/
int
ncmpi_mget_var_all(int           ncid, 
                   int           num, 
                   int           varids[],     /* [num] */
                   void         *bufs[],       /* [num] */
                   MPI_Offset    bufcounts[],  /* [num] */
                   MPI_Datatype  datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

/*----< ncmpi_mget_var1_all() >----------------------------------------------*/
int
ncmpi_mget_var1_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

/*----< ncmpi_mget_vara_all() >----------------------------------------------*/
int
ncmpi_mget_vara_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

/*----< ncmpi_mget_vars_all() >----------------------------------------------*/
int
ncmpi_mget_vars_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    MPI_Offset* const  strides[],    /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

/*----< ncmpi_mget_varm_all() >----------------------------------------------*/
int
ncmpi_mget_varm_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    MPI_Offset* const  strides[],    /* [num] */
                    MPI_Offset* const  imaps[],      /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

/*----< ncmpii_mgetput_varm() >-----------------------------------------------*/
static int
ncmpii_mgetput_varm(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[], /* [num] */
                    int                rw_flag,     /* WRITE_REQ or READ_REQ */
                    int                io_method)   /* COLL_IO or INDEP_IO */
{
    int i, status=NC_NOERR, *req_ids=NULL, *statuses=NULL;
    NC *ncp=NULL;

    CHECK_NCID
    if (rw_flag == WRITE_REQ)
        CHECK_WRITE_PERMISSION

    if (NC_indef(ncp)) return NC_EINDEFINE;

    /* check to see that the desired MPI file handle is opened */
    if (io_method == COLL_IO)
        CHECK_COLLECTIVE_FH
    else
        CHECK_INDEP_FH
  
    if (num > 0) {
        req_ids  = (int*) NCI_Malloc(2 * num * sizeof(int));
        statuses = req_ids + num;
    }

    /* for each request call ncmpi_igetput_varm() */
    for (i=0; i<num; i++) {
        NC_var *varp;
        MPI_Offset *start, *count;

        CHECK_VARID(varids[i], varp)

        if (starts == NULL) {         /* var */
            GET_FULL_DIMENSIONS
            status = ncmpii_igetput_varm(ncp, varp, start, count, NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag);
            if (varp->ndims > 0) NCI_Free(start);
        } else if (counts == NULL) {  /* var1 */
            GET_FULL_DIMENSIONS
            GET_ONE_COUNT
            status = ncmpii_igetput_varm(ncp, varp, starts[i], count, NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag);
            if (varp->ndims > 0) NCI_Free(count);
        } else if (strides == NULL) { /* vara */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i], NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag);
        } else if (imaps == NULL) {   /* vars */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag);
        } else {                      /* varm */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], imaps[i], bufs[i],
                                         bufcounts[i], datatypes[i],
                                         &req_ids[i], rw_flag);
        }
    }

    if (status != NC_NOERR)
        return status;

    if (io_method == COLL_IO)
        status = ncmpi_wait_all(ncid, num, req_ids, statuses);
    else
        status = ncmpi_wait(ncid, num, req_ids, statuses);
    if (status != NC_NOERR)
        return status;

    if (num > 0)
        NCI_Free(req_ids);

    for (i=0; i<num; i++)
        if (statuses[i] != NC_NOERR)
            return statuses[i];

    return NC_NOERR;
}
