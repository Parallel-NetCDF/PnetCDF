/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

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
    int i, status=NC_NOERR, *req_ids=NULL, *statuses=NULL, min_st;
    NC *ncp=NULL;

    /* check if ncid is valid */
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        /* must return the error now, parallel program might hang */
        return status;

    /* check if it is in define mode */
    if (NC_indef(ncp)) status = NC_EINDEFINE;

    /* check file write permission if this is write request */
    if (status == NC_NOERR) {
        if (rw_flag == WRITE_REQ && NC_readonly(ncp)) status = NC_EPERM;
    }
    /* check whether collective or independent mode */
    if (status == NC_NOERR) {
        if (io_method == INDEP_IO)
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->independent_fh),
                                        MPI_COMM_SELF, 0);
        else if (io_method == COLL_IO)
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->collective_fh),
                                        ncp->nciop->comm, 1);
        /* else if (io_method == INDEP_COLL_IO) */
    }
    if (io_method == COLL_IO)
        MPI_Allreduce(&status, &min_st, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
    else
        min_st = status;

    if (min_st != NC_NOERR)
        return status;
  
    if (num > 0) {
        req_ids  = (int*) NCI_Malloc(2 * num * sizeof(int));
        statuses = req_ids + num;
    }

    /* for each request call ncmpi_igetput_varm() */
    for (i=0; i<num; i++) {
        NC_var *varp;
        MPI_Offset *start, *count;

        varp = ncmpii_NC_lookupvar(ncp, varids[i]);
        if (varp == NULL) continue; /* invalid varid, skip this request */

        if (starts == NULL) {         /* var */
            GET_FULL_DIMENSIONS(start, count)
            status = ncmpii_igetput_varm(ncp, varp, start, count, NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag, 0);
            if (varp->ndims > 0) NCI_Free(start);
        } else if (counts == NULL) {  /* var1 */
            GET_FULL_DIMENSIONS(start, count)
            GET_ONE_COUNT(count)
            status = ncmpii_igetput_varm(ncp, varp, starts[i], count, NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag, 0);
            if (varp->ndims > 0) NCI_Free(count);
        } else if (strides == NULL) { /* vara */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i], NULL,
                                         NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag, 0);
        } else if (imaps == NULL) {   /* vars */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], NULL, bufs[i], bufcounts[i],
                                         datatypes[i], &req_ids[i], rw_flag, 0);
        } else {                      /* varm */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], imaps[i], bufs[i],
                                         bufcounts[i], datatypes[i],
                                         &req_ids[i], rw_flag, 0);
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
