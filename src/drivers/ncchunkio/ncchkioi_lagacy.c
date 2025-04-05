/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_get_var<kind>_all()        : dispatcher->get_var()
 * ncmpi_put_var<kind>_all()        : dispatcher->put_var()
 * ncmpi_get_var<kind>_<type>_all() : dispatcher->get_var()
 * ncmpi_put_var<kind>_<type>_all() : dispatcher->put_var()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncchkio_driver.h>
#include "ncchkio_internal.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

int
ncchkioi_get_var_old(NC_chk        *ncchkp,
              NC_chk_var       *varp,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int i, j, err;
    nc_type xtype;
    int *cstart, *cend, *ccord;
    int nb, bsize;
    int datavarid;
    int *bidx;
    int *tsize, *tssize, *tstart;
    int tpos;
    MPI_Datatype subarytype;
    char *rbuffer, *cbuffer;
    MPI_Offset cbsize;
    MPI_Offset **starts, **counts;

    // Boundary of chunks involved
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ccord = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        cstart[i] = start[i] / varp->chunkdim[i];
        if (stride == NULL){
            cend[i] = (start[i] + count[i] - 1) / varp->chunkdim[i];
        }
        else{
            cend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->chunkdim[i] + 1;
        }
    }
    
    // Number of chunks involved
    nb = 1;
    for(i = 0; i < varp->ndim; i++){
        nb *= cend[i] - cstart[i];
    }

    /* Use a varn call to read all compressed chunk involved 
     * Generate one request for each chunk
     */

    bidx = (int*)NCI_Malloc(sizeof(int) * nb);
    starts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nb);
    counts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nb);
    // Iterate through all chunks involved
    i = 0;
    cbsize = 0;
    memcpy(ccord, cstart, sizeof(int) * varp->ndim);
    for(i = 0; i < nb; i++){
        j = get_chunk_idx(varp, ccord);   
        bidx[i] = j; // chunk idx
        cbsize += varp->data_lens[j];  // total buffer size of compressed data
        starts[i] = varp->chunk_index + j;   // start of the chunk
        counts[i] = varp->data_lens + j; // count of the chunk

        // move on to next chunk
        ccord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (ccord[j] >= cend[j]){
                ccord[j - 1]++;
                ccord[j] = cstart[j];
            }
        }
    }

    // Allocate buffers
    cbuffer = (char*)NCI_Malloc(cbsize);  // Compressed data

    // Locate data var
    err = ncchkp->driver->get_var(ncchkp->ncp, varp->varid, NULL, NULL, NULL, NULL, &datavarid, 1, MPI_INT, reqMode); 
    if (err != NC_NOERR) return err;

    // read compressed data
    err = ncchkp->driver->get_varn(ncchkp->ncp, datavarid, nb, starts, counts, cbuffer, cbsize, MPI_BYTE, reqMode); 
    if (err != NC_NOERR) return err;

    // Decompression

    // Calculate chunk size
    // Original datatype
    err = ncchkp->driver->get_att(ncchkp->ncp, varp->varid, "_datatype", &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;

    // Calculate chunk size
    bsize = (int)NC_Type_size(xtype);
    for(i = 0; i < varp->ndim; i++){
        bsize *= varp->chunkdim[i];
    }

    // Allocate buffers
    rbuffer = NCI_Malloc(bsize * nb);  // Decompressed data

    // Decompress chunks
    cbsize = 0;
    for(i = 0; i < nb; i++){
        j = bidx[i];
        if (varp->data_lens[j] > 0){
            varp->filter_driver->decompress(cbuffer + cbsize, varp->data_lens[j], rbuffer + bsize * i, NULL, varp->ndim, varp->dimsize, ncmpii_nc2mpitype(xtype));
        }
        else{
            memset(rbuffer + bsize * i, 0, bsize);
        }
        cbsize += varp->data_lens[j];  // move to next chunk location
    }

    // Copy data into user buffer

    // Create datatype of querying domain in the decompressed domain
    tsize = NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = (cend[i] - cstart[i]) * varp->chunkdim[i];
        tssize[i] = (int)count[i];
        tstart[i] = start[i] % varp->chunkdim[i];
    }
    CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &subarytype);
    CHK_ERR_TYPE_COMMIT(&subarytype);

    // Pack data into user buffer
    tpos = 0;
    CHK_ERR_PACK(rbuffer, bsize * nb, subarytype, buf, bsize * nb, &tpos, ncchkp->comm);
    
    // Free datatype
    MPI_Type_free(&subarytype);

    return NC_NOERR;
}

int
ncchkioi_put_var_old(NC_chk        *ncchkp,
              NC_chk_var       *varp,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              void       *buf)
{
    int i, j, k, err;
    nc_type xtype;  // Variable data type in NC
    MPI_Datatype etype; // Variable element type in MPI
    int esize;  // Variable element size
    int *cstart, *cend, *ccord; // Bounding box for chunks overlapping my own write region
    int nb, bsize;  //number of chunks this process write to and chunk size
    int datavarid;  // Id of data variable
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int nmychunks, *mychunks;  // chunk count and id this process handles
    int *sendcounts, *sdispls;  // Send count and displacements in buffer
    int *recvcounts, *rdispls;  // Receive count and displacement in buffer
    int *packoff;   // Offset in mpi packing
    int *zipsize, *zdispls;  // Compressed count and displacement of my chunks in buffer
    int *zsize_local, *zsize_all;   // Compressed size of all chunks at local and global (all processes)
    int *zdispls_all;  // Compressed displacement of all chunks (all processes)
    int overlapsize;   // Size of overlapping region between a chunk and write region
    MPI_Datatype ptype;  // Pack datatype
    char *zbuf, *xbuf;  // Compressed and uncompressed data buffer
    char *sbuf, *rbuf;  // Send and receive buffer
    MPI_Offset **start_all, **count_all, **stride_all;  // Start, count, stride of all processes
    char name[128]; // Name of objects
    int zdimid;  // dimension id for compressed data variable
    MPI_Offset **zstarts, **zcounts;    // Starts and counts in the varn call for compressed data

    // Original datatype and size
    err = ncchkp->driver->get_att(ncchkp->ncp, varp->varid, "_datatype", &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;
    esize = NC_Type_size(xtype);
    etype = ncmpii_nc2mpitype(xtype);

    // Calculate chunk size
    bsize = esize;
    for(i = 0; i < varp->ndim; i++){
        bsize *= varp->chunkdim[i];
    }

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    /*
     * Gather start, count, stride to all processes
     */

    // Allocate buffer

    start_all = NCI_Malloc(sizeof(MPI_Offset*) * ncchkp->np);
    count_all = NCI_Malloc(sizeof(MPI_Offset*) * ncchkp->np);
    stride_all = NCI_Malloc(sizeof(MPI_Offset*) * ncchkp->np);

    start_all[0] = NCI_Malloc(sizeof(MPI_Offset) * ncchkp->np * varp->ndim);
    count_all[0] = NCI_Malloc(sizeof(MPI_Offset) * ncchkp->np * varp->ndim);
    stride_all[0] = NCI_Malloc(sizeof(MPI_Offset) * ncchkp->np * varp->ndim);

    for(i = 1; i < ncchkp->np; i++){
        start_all[i] = start_all[0] + i * varp->ndim;
        count_all[i] = count_all[0] + i * varp->ndim;
        stride_all[i] = stride_all[0] + i * varp->ndim;
    }

    // Call allgather

    err = MPI_Allgather(start, varp->ndim, MPI_LONG_LONG_INT, start_all[0], varp->ndim, MPI_LONG_LONG_INT, ncchkp->comm);
    if (err != MPI_SUCCESS){
        err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
        DEBUG_RETURN_ERROR(err);
    }

    if (count != NULL){
        err = MPI_Allgather(count, varp->ndim, MPI_LONG_LONG_INT, count_all[0], varp->ndim, MPI_LONG_LONG_INT, ncchkp->comm);
        if (err != MPI_SUCCESS){
            err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
            DEBUG_RETURN_ERROR(err);
        }
    }

    if (stride != NULL){
        err = MPI_Allgather(stride, varp->ndim, MPI_LONG_LONG_INT, stride_all[0], varp->ndim, MPI_LONG_LONG_INT, ncchkp->comm);
        if (err != MPI_SUCCESS){
            err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
            DEBUG_RETURN_ERROR(err);
        }
    }

    /* 
     * Now, we need to send data to the chunk owner as well as receive data for our own chunk
     */

    // First, compute chunk boundary, find overlapping chunks
    cstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    ccord = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    cend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        cstart[i] = start[i] / varp->chunkdim[i];
        if (stride == NULL){
            cend[i] = (start[i] + count[i] - 1) / varp->chunkdim[i] + 1;
        }
        else{
            cend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->chunkdim[i] + 1;
        }
    }

    // Calculate the amount we need to send to other process
    sendcounts = (int*)NCI_Malloc(sizeof(int) * ncchkp->np);
    sdispls = (int*)NCI_Malloc(sizeof(int) * ncchkp->np);
    packoff = (int*)NCI_Malloc(sizeof(int) * ncchkp->np);
    memset(sendcounts, 0, sizeof(int) * ncchkp->np);
    memset(packoff, 0, sizeof(int) * ncchkp->np);

    // Iterate through all chunks involved to count send size
    i = 0;
    overlapsize = 0;
    memcpy(ccord, cstart, sizeof(int) * varp->ndim);
    while(ccord[0] < cend[0]){
        j = varp->chunk_owner[get_chunk_idx(varp, ccord)];    
        
        // Overlapping size of this chunk
        overlapsize = get_chunk_overlap(varp, ccord, start, count, stride, tstart, tssize);
        sendcounts[j] += overlapsize;

        // move on to next chunk
        ccord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (ccord[j] >= cend[j]){
                ccord[j - 1]++;
                ccord[j] = cstart[j];
            }
        }
    }

    // Buffer displacement
    sdispls[0] = 0;
    for(i = 1; i < ncchkp->np; i++){
        sdispls[i] = sendcounts[i - 1] + sdispls[i - 1];
    }

    // Allocate send buffer
    sbuf = (char*)NCI_Malloc(sdispls[ncchkp->np - 1] + sendcounts[ncchkp->np - 1]);

    // Pack data into send buffer
    
    // Iterate through all chunks involved again, this time actually pack the data
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = (int)count[i];
    }
    i = 0;
    overlapsize = 0;
    memcpy(ccord, cstart, sizeof(int) * varp->ndim);
    while(ccord[0] < cend[0]){
        j = varp->chunk_owner[get_chunk_idx(varp, ccord)];   
        
        // Overlapping region of this chunk
        get_chunk_overlap(varp, ccord, start, count, stride, tstart, tssize);
        for(k = 0; k < varp->ndim; k++){
            tstart[k] -= (int)start[k];
        }

        // Pack type
        CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &ptype);
        CHK_ERR_TYPE_COMMIT(&ptype);
        
        // Pack data
        CHK_ERR_PACK(buf, 1, ptype, sbuf + sdispls[j], sendcounts[j], packoff + j, ncchkp->comm);
    
        // Free datatype
        MPI_Type_free(&ptype);

        // move on to next chunk
        ccord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (ccord[j] >= cend[j]){
                ccord[j - 1]++;
                ccord[j] = cstart[j];
            }
        }
    }

    /*
     * Determine chunk ownership
     * Find my chunks
     */
    nmychunks = 0;
    for(i = 0; i < varp->nchunk; i++){
        if (varp->chunk_owner[i] == ncchkp->rank){
            nmychunks++;
        }
    }

    // Gather chunk id this process handled to prevent a search in the future
    mychunks = (int*)NCI_Malloc(sizeof(int) * nmychunks);
    nmychunks = 0;
    for(i = 0; i < varp->nchunk; i++){
        if (varp->chunk_owner[i] == ncchkp->rank){
            mychunks[nmychunks] = i;
            nmychunks++;
        }
    }

    /* 
     * Compute size to receive
     * We only need size here, packing will happen after receving
     */

    // Calculate the amount we need to receive from other process
    recvcounts = (int*)NCI_Malloc(sizeof(int) * ncchkp->np);
    rdispls = (int*)NCI_Malloc(sizeof(int) * ncchkp->np);
    memset(recvcounts, 0, sizeof(int) * ncchkp->np);
    memset(packoff, 0, sizeof(int) * ncchkp->np);
    for(i = 0; i < varp->nchunk; i++){
        if (varp->chunk_owner[i] == ncchkp->rank){
            get_chunk_cord(varp, i, ccord);

            for(j = 0; j < ncchkp->np; j++){
                // Overlapping region of this chunk
                get_chunk_overlap(varp, ccord, start_all[j], count_all[j], stride_all[j], tstart, tssize);

                overlapsize = esize;
                for(k = 0; k < varp->ndim; k++){
                    overlapsize *= tssize[k];
                }
                recvcounts[j] += overlapsize;
            }
        }
    }

    // Buffer displacement
    rdispls[0] = 0;
    for(i = 1; i < ncchkp->np; i++){
        rdispls[i] = recvcounts[i - 1] + rdispls[i - 1];
    }

    // Allocate receive buffer
    rbuf = (char*)NCI_Malloc(rdispls[ncchkp->np - 1] + recvcounts[ncchkp->np - 1]);

    // Send the data to destination
    MPI_Alltoallv(sbuf, sendcounts, sdispls, MPI_BYTE, rbuf, recvcounts, rdispls, MPI_BYTE, ncchkp->comm);

/*
#ifdef PNETCDF_DEBUG
    if (ncchkp->rank == 0){
        printf("Rank %d: sendcount = {", ncchkp->rank);
        for(i = 0; i < ncchkp->np; i++){
            printf("%d, ", sendcounts[i]);
        }
        printf("}, sdispls = {");
        for(i = 0; i < ncchkp->np; i++){
            printf("%d, ", sdispls[i]);
        }
        printf("}, recvcounts = {");
        for(i = 0; i < ncchkp->np; i++){
            printf("%d, ", recvcounts[i]);
        }
        printf("}, rdispls = {");
        for(i = 0; i < ncchkp->np; i++){
            printf("%d, ", rdispls[i]);
        }
        printf("}, sbuf = {");
        for(i = 0; i < sdispls[ncchkp->np - 1] + sendcounts[ncchkp->np - 1]; i++){
            printf("%x ", sbuf[i]);
        }
        printf("}, rbuf = {");
        for(i = 0; i < rdispls[ncchkp->np - 1] + recvcounts[ncchkp->np - 1]; i++){
            printf("%x ", rbuf[i]);
        }
        printf("}\n");
        fflush(stdout);
    }
#endif
*/

    /*
     * Next step is to pack data to chunk buffer
     */

    // Allocate buffer
    xbuf = (char*)NCI_Malloc(nmychunks * bsize);

    // Main array is the whole chunk
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->chunkdim[i];
    }

    // Pack data
    memset(packoff, 0, sizeof(int) * ncchkp->np);
    for(i = 0; i < nmychunks; i++){
        get_chunk_cord(varp, mychunks[i], ccord);

        for(j = 0; j < ncchkp->np; j++){
            // Overlapping region of this chunk
            overlapsize = get_chunk_overlap(varp, ccord, start_all[j], count_all[j], stride_all[j], tstart, tssize);

            if (overlapsize > 0){
                // Overlap size
                //overlapsize = esize;
                //for(k = 0; k < varp->ndim; k++){
                //    overlapsize *= tssize[k];
                //}

                // The chunk is the main array, overlapping region is the subarray
                for(k = 0; k < varp->ndim; k++){
                    tstart[k] -= ccord[k] * varp->chunkdim[k];
                }

                // Pack type
                CHK_ERR_TYPE_CREATE_SUBARRAY(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &ptype);
                CHK_ERR_TYPE_COMMIT(&ptype);

                // Pack data
                CHK_ERR_UNPACK(rbuf + rdispls[j], overlapsize, packoff + j, xbuf + bsize * i, 1, ptype, ncchkp->comm);

                // Free datatype
                MPI_Type_free(&ptype);
            }
        }
    }

/*
#ifdef PNETCDF_DEBUG
    if (ncchkp->rank == 0){
        printf("Rank %d: xbuf = {", ncchkp->rank);
        for(i = 0; i < nmychunks * bsize; i++){
            printf("%x ", xbuf[i]);
        }
        printf("}\n");
        fflush(stdout);
    }
#endif
*/

    /* 
     * The buffer is now filled with data coming from all processes, it's time to compress
     */

    // compressed size and displacement
    zipsize = (int*)NCI_Malloc(sizeof(int) * nmychunks);
    zdispls = (int*)NCI_Malloc(sizeof(int) * (nmychunks + 1));
    memset(zipsize, 0, sizeof(int) * nmychunks);
    memset(zdispls, 0, sizeof(int) * (nmychunks + 1));

    // Calculate compressed data size
    for(i = 0; i < nmychunks; i++){
        // Calculate compressed size
        // This is just estimate
        varp->filter_driver->compress(xbuf + bsize * i, bsize, NULL, zipsize + i, varp->ndim, varp->chunkdim, etype);
    }

    // Calculate total size
    for(i = 1; i < nmychunks; i++){
        zdispls[0] += zipsize[i];
    }

    // Allocate buffer
    zbuf = (char*)NCI_Malloc(zdispls[0]);

    // Perform real compression
    for(i = 0; i < nmychunks; i++){
        // Compressed the data
        // We get real size here
        varp->filter_driver->compress(xbuf + bsize * i, bsize, zbuf + zdispls[i], zipsize + i, varp->ndim, varp->chunkdim, etype);
        
        // Calculate offset
        zdispls[i + 1] = zdispls[i] + zipsize[i];
    }

/*
#ifdef PNETCDF_DEBUG
    if (ncchkp->rank == 0){
        printf("Rank %d: zipsize = {", ncchkp->rank);
        for(i = 0; i < nmychunks; i++){
            printf("%x ", zipsize[i]);
        }
        printf("}, zdispls = {");
        for(i = 0; i < nmychunks; i++){
            printf("%d, ", zdispls[i]);
        }
        printf("}, zbuf = {");
        for(i = 0; i < zdispls[nmychunks- 1] + zipsize[nmychunks - 1]; i++){
            printf("%x ", zbuf[i]);
        }
        printf("}\n");
        fflush(stdout);
    }
#endif
*/

    /*
     * Now it is time for a collective write
     * We start by syncing compressed size on all processes
     * Then, we can create variable large enough to store compressed data
     * Finally, we do collective write to store the data
     */

    // First sync on compressed chunk size
    // We use a all MAX reduce on all chunks
    // An alternative is to allgather and unpack the info

    // Allocate buffer
    zsize_local = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    zsize_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    zdispls_all = (int*)NCI_Malloc(sizeof(int) * varp->nchunk);
    memset(zsize_local, 0, sizeof(int) * varp->nchunk);
    memset(zsize_all, 0, sizeof(int) * varp->nchunk);
    memset(zdispls_all, 0, sizeof(int) * varp->nchunk);

    // Fill up local size
    for(i = 0; i < nmychunks; i++){
        zsize_local[mychunks[i]] = zipsize[i];
    }

    // All reduce
    CHK_ERR_ALLREDUCE(zsize_local, zsize_all, varp->nchunk, MPI_INT, MPI_MAX, ncchkp->comm);

    // Calculate variable displacement
    zdispls_all[0] = 0;
    for(i = 1; i < varp->nchunk; i++){
        zdispls_all[i] = zsize_all[i - 1] + zdispls_all[i - 1];
    }

/*
#ifdef PNETCDF_DEBUG
    if (ncchkp->rank == 0){
        printf("Rank %d: zsize_all = {", ncchkp->rank);
        for(i = 0; i < varp->nchunk; i++){
            printf("%x ", zsize_all[i]);
        }
        printf("}, zdispls_all = {");
        for(i = 0; i < varp->nchunk; i++){
            printf("%d, ", zdispls_all[i]);
        }
        printf("}, varid = { %d", varp->varid);
        printf("}, datavarid = { %d", varp->datavarid);
        printf("}\n");
        fflush(stdout);
    }
#endif
*/

    // Enter redefine mode
    ncchkp->driver->redef(ncchkp->ncp);

    // Define dimension  for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid);
    err = ncchkp->driver->def_dim(ncchkp->ncp, name, zdispls_all[varp->nchunk - 1] + zsize_all[varp->nchunk - 1], &zdimid);
    if (err != NC_NOERR) return err;

    // Define variable
    sprintf(name, "_compressed_data_%d", varp->varid);
    err = ncchkp->driver->def_var(ncchkp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));
    if (err != NC_NOERR) return err;

    // Record offset in data variable
    err = ncchkp->driver->put_att(ncchkp->ncp, varp->varid, "_chunkoffset", NC_INT, varp->nchunk, zdispls_all, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Switch to data mode
    err = ncchkp->driver->enddef(ncchkp->ncp);
    if (err != NC_NOERR) return err;

    //Now, we generate a varn call to write out compressed data
    zstarts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nmychunks);
    zcounts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * nmychunks);
    zstarts[0] = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * nmychunks);
    zcounts[0] = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * nmychunks);
    for(i = 0; i < nmychunks; i++){
        zstarts[i] = zstarts[0] + i;
        zcounts[i] = zcounts[0] + i;
        zstarts[i][0] = zdispls_all[mychunks[i]];
        zcounts[i][0] = zsize_all[mychunks[i]];
    }
    err = ncchkp->driver->put_varn(ncchkp->ncp, varp->datavarid, nmychunks, zstarts, zcounts, zbuf, zdispls[nmychunks - 1] + zipsize[nmychunks - 1], MPI_UNSIGNED_CHAR, NC_REQ_WR | NC_REQ_BLK | NC_REQ_FLEX | NC_REQ_COLL);
    if (err != NC_NOERR) return err;

    // Record datavar id
    err = ncchkp->driver->put_var(ncchkp->ncp, varp->varid, NULL, NULL, NULL, NULL, &(varp->datavarid), 1, MPI_INT, NC_REQ_WR | NC_REQ_BLK | NC_REQ_FLEX | NC_REQ_COLL);
    if (err != NC_NOERR) return err;
    
    //  Free up buffers
    NCI_Free(cstart);
    NCI_Free(cend);
    NCI_Free(ccord);
    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);
    NCI_Free(mychunks);
    NCI_Free(sendcounts);
    NCI_Free(sdispls);
    NCI_Free(recvcounts);
    NCI_Free(rdispls);
    NCI_Free(packoff);
    NCI_Free(zipsize);
    NCI_Free(zdispls);
    NCI_Free(zsize_local);
    NCI_Free(zsize_all);
    NCI_Free(zdispls_all);
    NCI_Free(zbuf);
    NCI_Free(xbuf);
    NCI_Free(sbuf);
    NCI_Free(rbuf);
    NCI_Free(start_all[0]);
    NCI_Free(count_all[0]);
    NCI_Free(stride_all[0]);
    NCI_Free(start_all);
    NCI_Free(count_all);
    NCI_Free(stride_all);
    NCI_Free(zstarts[0]);
    NCI_Free(zcounts[0]);
    NCI_Free(zstarts);
    NCI_Free(zcounts);

    return NC_NOERR;
}

void profile(){
        /* Profiling information */
    ncchkp->profile.total_data += t9 - t0;
    ncchkp->profile.total_meta += t9 - t0;
    ncchkp->profile.max_buffer += t9 - t0;
    ncchkp->profile.total_time += t9 - t0;
    ncchkp->profile.cb_time += t9 - t0;
    ncchkp->profile.io_time += t9 - t0;
    
    ncchkp->profile.cb_init_time += t9 - t0;    // Calculate number of req
    ncchkp->profile.cb_sync_time += t9 - t0;    // Syncing number of req
    ncchkp->profile.cb_pack_req_time += t9 - t0;    // Pack request and reply
    ncchkp->profile.cb_pack_rep_time += t9 - t0;    // Pack request and reply
    ncchkp->profile.cb_unpack_req_time += t9 - t0;    // Unpack incoming request
    ncchkp->profile.cb_unpack_rep_time += t9 - t0;    // Unpack incoming request
    ncchkp->profile.cb_send_req_time += t9 - t0;    // Posting and waiting send
    ncchkp->profile.cb_send_rep_time += t9 - t0;    // Posting and waiting send
    ncchkp->profile.cb_recv_req_time += t9 - t0;    // Time posting and waiting recv
    ncchkp->profile.cb_recv_rep_time += t9 - t0;    // Time posting and waiting recv
    ncchkp->profile.cb_self_time += t9 - t0;    // Time handling our own data

    ncchkp->profile.io_wr_time += t9 - t0;
    ncchkp->profile.io_rd_time += t9 - t0;
    ncchkp->profile.io_com_time += t9 - t0;
    ncchkp->profile.io_decom_time += t9 - t0;
    ncchkp->profile.io_sync_time += t9 - t0;
}