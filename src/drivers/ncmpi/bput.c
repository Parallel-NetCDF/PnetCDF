/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "ncmpidtype.h"
#include "macro.h"

/*----< ncmpii_buffer_attach() >---------------------------------------------*/
int
ncmpii_buffer_attach(void       *ncdp,
                     MPI_Offset  bufsize)
{
    NC *ncp=(NC*)ncdp;

    if (bufsize <= 0) DEBUG_RETURN_ERROR(NC_ENULLBUF)

    /* check if the buffer has been previously attached
     * note that in nc.c, the NC object is allocated with calloc, so
     * abuf should be initialized to NULL then
     */
    if (ncp->abuf != NULL) DEBUG_RETURN_ERROR(NC_EPREVATTACHBUF)

    ncp->abuf = (NC_buf*) NCI_Malloc(sizeof(NC_buf));

    ncp->abuf->size_allocated = bufsize;
    ncp->abuf->size_used = 0;
    ncp->abuf->table_size = NC_ABUF_DEFAULT_TABLE_SIZE;
    ncp->abuf->occupy_table = (NC_buf_status*)
               NCI_Calloc(NC_ABUF_DEFAULT_TABLE_SIZE, sizeof(NC_buf_status));
    ncp->abuf->tail = 0;
    ncp->abuf->buf = NCI_Malloc((size_t)bufsize);
    return NC_NOERR;
}

/*----< ncmpii_buffer_detach() >---------------------------------------------*/
int
ncmpii_buffer_detach(void *ncdp)
{
    int  i;
    NC  *ncp=(NC*)ncdp;

    /* check if the buffer has been previously attached */
    if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)

    /* this API assumes users are responsible for no pending bput */
    for (i=0; i<ncp->numPutReqs; i++) {
        if (ncp->put_list[i].abuf_index >= 0) /* check for a pending bput */
            DEBUG_RETURN_ERROR(NC_EPENDINGBPUT)
            /* return now, so users can call wait and try detach again */
    }

    NCI_Free(ncp->abuf->buf);
    NCI_Free(ncp->abuf->occupy_table);
    NCI_Free(ncp->abuf);
    ncp->abuf = NULL;

    return NC_NOERR;
}

#ifdef THIS_SEEMS_OVER_DONE_IT
/*----< ncmpii_buffer_detach() >---------------------------------------------*/
/* mimic MPI_Buffer_detach()
 * Note from MPI: Even though the 'bufferptr' argument is declared as
 * 'void *', it is really the address of a void pointer.
 */
int
ncmpii_buffer_detach(void       *ncdp,
                     void       *bufptr,
                     MPI_Offset *bufsize)
{
    int  i;
    NC  *ncp=(NC*)ncdp;

    /* check if the buffer has been previously attached */
    if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)

    /* check MPICH2 src/mpi/pt2pt/bsendutil.c for why the bufptr is void* */
    *(void **)bufptr = ncp->abuf->buf;
    *bufsize         = ncp->abuf->size_allocated;

    /* this API assumes users are responsible for no pending bput when called */
    for (i=0; i<ncp->numPutReqs; i++) {
        if (ncp->put_list[i].abuf_index >= 0) /* check for a pending bput */
            DEBUG_RETURN_ERROR(NC_EPENDINGBPUT)
            /* return now, so users can call wait and try detach again */
    }

    NCI_Free(ncp->abuf->occupy_table);
    NCI_Free(ncp->abuf);
    ncp->abuf = NULL;

    return NC_NOERR;
}
#endif

/*----< ncmpii_bput_var() >--------------------------------------------------*/
int
ncmpii_bput_var(void             *ncdp,
                int               varid,
                const MPI_Offset *start,
                const MPI_Offset *count,
                const MPI_Offset *stride,
                const MPI_Offset *imap,
                const void       *buf,
                MPI_Offset        bufcount,
                MPI_Datatype      buftype,
                int              *reqid,
                int               api_kind,
                nc_type           itype)
{
    int         err;
    NC         *ncp=(NC*)ncdp;
    NC_var     *varp=NULL;
    MPI_Offset *_start, *_count;

    if (reqid != NULL) *reqid = NC_REQ_NULL;

    err = ncmpii_sanity_check(ncp, varid, start, count, stride,
                              bufcount, buftype, api_kind, (itype=NC_NAT),
                              0, WRITE_REQ, NONBLOCKING_IO, &varp);
    if (err != NC_NOERR) return err;

    if (ncp->abuf == NULL) DEBUG_RETURN_ERROR(NC_ENULLABUF)

    _start = (MPI_Offset*)start;
    _count = (MPI_Offset*)count;
         if (api_kind == API_VAR)  GET_FULL_DIMENSIONS(_start, _count)
    else if (api_kind == API_VAR1) GET_ONE_COUNT(_count)

    err = ncmpii_igetput_varm(ncp, varp, _start, _count, stride, imap,
                              (void*)buf, bufcount, buftype,
                              reqid, WRITE_REQ, 1, 0);

         if (api_kind == API_VAR)  NCI_Free(_start);
    else if (api_kind == API_VAR1) NCI_Free(_count);

    return err;
}
