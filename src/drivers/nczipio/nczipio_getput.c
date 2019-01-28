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

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nczipio_driver.h>

static int get_block_idx(NC_var *varp, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->stripesize[i - 1] + cord[i];
    }

    return ret;
}

static int get_block_cord(NC_var *varp, int idx, int* cord){
    int i, ret;
    
    ret = cord[0];
    for(i = 1; i < varp->ndim; i++){
        ret = ret * varp->stripesize[i - 1] + cord[i];
    }

    for(i = varp->ndim - 1; i > 0; i--){
        cord[i] = idx % varp->stripesize[i - 1];
        idx /= varp->stripesize[i - 1];
    }
    cord[0] = idx;

    return 0;
}

static int get_block_overlap(NC_var *varp, int* cord, MPI_Offset *start, MPI_Offset *count, MPI_Offset *stride, int *ostart, int *ocount){
    int i, ret;
    
    for(i = 0; i < varp->ndim; i++){
        ostart[i] = max(start[i], cord[i] * varp->stripesize[i]);
        ocount[i] = min(start[i] + count[i], (cord[i] + 1) * varp->stripesize[i]) - ostart[i];
        if (ocount[i] < 0){
            ocount[i] = 0;
        }
    }

    return 0;
}

int
nczipio_get_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int i, err;
    NC_zip *nczipp = (NC_zip*)ncdp;
    NC_var *varp;
    nc_type xtype;
    int *bstart, *bend, *bcord;
    int nb, bsize;
    int datavarid;
    int *tsize, *tssize, *tstart;
    int tpos;
    MPI_Datatype subarytype;
    char *rbuffer, *cbuffer;
    MPI_Offset cbsize;
    MPI_Offset **starts, **counts;

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    // Boundary of blocks involved
    *bstart = NCI_Malloc(sizeof(int) * varp->ndim);
    *bcord = NCI_Malloc(sizeof(int) * varp->ndim);
    *bend = NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        bstart[i] = start[i] / varp->stripesize;
        if (stride == NULL){
            bend[i] = (start[i] + count[i] - 1) / varp->stripesize;
        }
        else{
            bend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->stripesize + 1;
        }
    }
    
    // Number of blocks involved
    nb = 1;
    for(i = 0; i < varp->ndim; i++){
        nb *= bend[i] - bstart[i];
    }

    /* Use a varn call to read all compressed block involved 
     * Generate one request for each block
     */

    *bidx = NCI_Malloc(sizeof(int) * nb);
    **starts = NCI_Malloc(sizeof(MPI_Offset*) * nb);
    **counts = NCI_Malloc(sizeof(MPI_Offset*) * nb);
    // Iterate through all blocks involved
    i = 0;
    cbsize = 0;
    memcpy(bcord, bstart, sizeof(int) * varp->ndim);
    for(i = 0; i < nb; i++){
        j = get_block_idx(varp, bcord);   
        bidx[i] = j; // block idx
        cbsize += varp->lens[j];  // total buffer size of compressed data
        starts[i] = varp->offset + j;   // start of the block
        counts[i] = varp->lens + j; // count of the block

        // move on to next block
        bcord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (bcord[j] >= bend[j]){
                bcord[j - 1]++;
                bcord[j] = bstart[j];
            }
        }
    }

    // Allocate buffers
    *cbuffer = NCI_Malloc(cbsize);  // Compressed data

    // Locate data var
    err = nczipp->driver->get_var(nczipp->ncp, varp->varid, NULL, NULL, NULL, NULL, &datavarid, 1, MPI_INT, reqMode); 
    if (err != NC_NOERR) return err;

    // read compressed data
    err = nczipp->driver->get_varn(nczipp->ncp, datavarid, nb, starts, counts, cbuffer, cbsize, MPI_BYTE, reqMode); 
    if (err != NC_NOERR) return err;

    // Decompression

    // Calculate block size
    // Original datatype
    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;

    // Calculate block size
    bsize = (int)NC_Type_size(xtype);
    for(i = 0; i < nblocks; i++){
        bsize *= varp->stripesize[i];
    }

    // Allocate buffers
    *rbuffer = NCI_Malloc(bsize * nb);  // Decompressed data

    // Decompress blocks
    cbsize = 0;
    for(i = 0; i < nb; i++){
        j = bidx[i];
        if (varp->lens[j] > 0){
            nczipp->zip->decompress(cbuffer + cbsize, varp->lens[j], rbuffer + bsize * i, NULL, carp->ndim, varp->dimsize, ncmpii_nc2mpitype(xtype));
        }
        else{
            memset(rbuffer + bsize * i, 0, bsize);
        }
        cbsize += varp->lens[j];  // move to next block location
    }

    // Copy data into user buffer

    // Create datatype of querying domain in the decompressed domain
    *tsize = NCI_Malloc(sizeof(int) * varp->ndim);
    *tssize = NCI_Malloc(sizeof(int) * varp->ndim);
    *tstart = NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = (bend[i] - bstart[i]) * varp->stripesize[i];
        tssize[i] = (int)count[i];
        tstart[i] = start[i] % varp->stripesize[i];
    }
    MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &subarytype);

    // Pack data into user buffer
    tpos = 0;
    MPI_Pack(rbuffer, bsize * nb, subarytype, buf, bsize * nb, &tpos, nczipp->comm);

    return NC_NOERR;
}

int
nczipio_get_varn(void              *ncdp,
               int                varid,
               int                num,
               MPI_Offset* const *starts,
               MPI_Offset* const *counts,
               void              *buf,
               MPI_Offset         bufcount,
               MPI_Datatype       buftype,
               int                reqMode)
{

    // Original datatype
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;

    err = nczipp->driver->get_var(nczipp->ncp, varid, start, count, stride, imap,
                               buf, bufcount, buftype, reqMode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
nczipioi_put_var(NC_zip             *nczipp,
              NC_var               *varp,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const void       *buf)
{
    int i, err;
    nc_type xtype;  // Variable data type in NC
    MPI_Datatype etype; // Variable element type in MPI
    int esize;  // Variable element size
    int *bstart, *bend, *bcord; // Bounding box for blocks overlapping my own write region
    int nb, bsize;  //number of blocks this process write to and block size
    int datavarid;  // Id of data variable
    int *tsize, *tssize, *tstart;   // Size for sub-array type
    int nmyblocks, *myblocks;  // Block count and id this process handles
    int *sendcounts, *sdispls;  // Send count and displacements in buffer
    int *recvcounts, *rdispls;  // Receive count and displacement in buffer
    int *zipsize, *zdispls;  // Compressed count and displacement of my blocks in buffer
    int *zsize_local, *zsize_all;   // Compressed size of all blocks at local and global (all processes)
    int *zdispls_all;  // Compressed displacement of all blocks (all processes)
    int overlapsize;   // Size of overlapping region between a block and write region
    MPI_Datatype ptype;  // Pack datatype
    char *zbuf, *xbuf;  // Compressed and uncompressed data buffer
    char *sbuf, *rbuf;  // Send and receive buffer
    MPI_Offset **start_all, **count_all, **stride_all;  // Start, count, stride of all processes
    char name[128]; // Name of objects
    int zdimid;  // dimension id for compressed data variable
    MPI_Offset **zstarts, **zcounts;    // Starts and counts in the varn call for compressed data

    if (varid < 0 || varid >= nczipp->vars.cnt){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }
    varp = nczipp->vars.data + varid;

    // Original datatype and size
    err = nczipp->driver->get_att(nczipp->ncp, varp->varid, "_datatype", NC_INT, 1, &xtype, MPI_INT); 
    if (err != NC_NOERR) return err;
    esize = NC_Type_size(xtype);
    etype = ncmpii_nc2mpitype(xtype);

    // Calculate block size
    bsize = esize;
    for(i = 0; i < varp->ndim; i++){
        bsize *= varp->stripesize[i];
    }

    // Allocate buffering for overlaping index
    tsize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tssize = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    tstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);

    /*
     * Gather start, count, stride to all processes
     */

    // Allocate buffer

    start_all = NCI_Malloc(sizeof(MPI_Offset*) * nczipp->np);
    count_all = NCI_Malloc(sizeof(MPI_Offset*) * nczipp->np);
    stride_all = NCI_Malloc(sizeof(MPI_Offset*) * nczipp->np);

    start_all[0] = NCI_Malloc(sizeof(MPI_Offset) * nczipp->np * varp->ndim);
    count_all[0] = NCI_Malloc(sizeof(MPI_Offset) * nczipp->np * varp->ndim);
    stride_all[0] = NCI_Malloc(sizeof(MPI_Offset) * nczipp->np * varp->ndim);

    for(i = 1; i < nczipp->np; i++){
        start_all[i] = start_all[0] + i * varp->ndim;
        count_all[i] = count_all[0] + i * varp->ndim;
        stride_all[i] = stride_all[0] + i * varp->ndim;
    }

    // Call allgather

    err = MPI_Allgather(start, varp->ndim, MPI_LONG_LONG_INT, start_all[0], nczipp->np * varp->ndim, MPI_LONG_LONG_INT, nczipp->comm);
    if (err != MPI_SUCCESS){
        err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
        DEBUG_RETURN_ERROR(err);
    }

    err = MPI_Allgather(count, varp->ndim, MPI_LONG_LONG_INT, count_all[0], nczipp->np * varp->ndim, MPI_LONG_LONG_INT, nczipp->comm);
    if (err != MPI_SUCCESS){
        err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
        DEBUG_RETURN_ERROR(err);
    }

    err = MPI_Allgather(stride, varp->ndim, MPI_LONG_LONG_INT, stride_all[0], nczipp->np * varp->ndim, MPI_LONG_LONG_INT, nczipp->comm);
    if (err != MPI_SUCCESS){
        err = ncmpii_error_mpi2nc(err, "MPI_Allgather");
        DEBUG_RETURN_ERROR(err);
    }

    /* 
     * Now, we need to send data to the block owner as well as receive data for our own block
     */

    // First, compute block boundary, find overlapping blocks
    bstart = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    bcord = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    bend = (int*)NCI_Malloc(sizeof(int) * varp->ndim);
    for(i = 0; i < varp->ndim; i++){
        bstart[i] = start[i] / varp->stripesize;
        if (stride == NULL){
            bend[i] = (start[i] + count[i] - 1) / varp->stripesize;
        }
        else{
            bend[i] = (start[i] + (count[i] - 1) * stride[i]) / varp->stripesize + 1;
        }
    }

    // Calculate the amount we need to send to other process
    sendcounts = (int*)NCI_Malloc(sizeof(int) * nczipp->np);
    sdispls = (int*)NCI_Malloc(sizeof(int) * nczipp->np);
    packoff = (int*)NCI_Malloc(sizeof(int) * nczipp->np);
    memset(sendcounts, 0, sizeof(int) * nczipp->np);
    memset(packoff, 0, sizeof(int) * nczipp->np);

    // Iterate through all blocks involved to count send size
    i = 0;
    overlapsize = 0;
    memcpy(bcord, bstart, sizeof(int) * varp->ndim);
    while(bcord[0] < bend[0]){
        j = get_block_idx(varp, bcord);   
        
        // Overlapping size of this block
        get_block_overlap(varp, bcord, start, count, stride, tstart, tssize);
        overlapsize = esize;
        for(k = 0; k < varp->ndim; k++){
            overlapsize *= tssize[k];
        }
        sendcounts[j] += overlapsize;

        // move on to next block
        bcord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (bcord[j] >= bend[j]){
                bcord[j - 1]++;
                bcord[j] = bstart[j];
            }
        }
    }

    // Buffer displacement
    sdispls[0] = 0;
    for(i = 1; i < nczipp->np; i++){
        sdispls[i] = sendcounts[i - 1] + sdispls[i - 1];
    }

    // Allocate send buffer
    sbuf = (char*)NCI_Malloc(sdispls[ncaipp->np - 1] + sendcounts[ncaipp->np - 1]);

    // Pack data into send buffer
    
    // Iterate through all blocks involved again, this time actually pack the data
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = (int)count[i];
    }
    i = 0;
    overlapsize = 0;
    memcpy(bcord, bstart, sizeof(int) * varp->ndim);
    while(bcord[0] < bend[0]){
        j = get_block_idx(varp, bcord);   
        
        // Overlapping region of this block
        get_block_overlap(varp, bcord, start, count, stride, tstart, tssize);
        for(k = 0; k < varp->ndim; k++){
            tstart[k] -= (int)start[k];
        }

        // Pack type
        MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &ptype);
        
        // Pack data
        MPI_Pack(buf, 1, ptype, sbuf + sdispls[j], sendcounts[j], packoff + j, nczipp->comm);

        // move on to next block
        bcord[varp->ndim - 1]++;
        for(j = varp->ndim - 1; j > 0; j--){
            if (bcord[j] >= bend[j]){
                bcord[j - 1]++;
                bcord[j] = bstart[j];
            }
        }
    }

    /*
     * Determine block ownership
     * Find my blocks
     */
    nmyblocks = 0;
    for(i = 0; i < varp->nblocks; i++){
        if (varp->owner[i] == nczipp->rank){
            nmyblocks++;
        }
    }

    // Gather block id this process handled to prevent a search in the future
    myblocks = (int*)NCI_Malloc(sizeof(int) * nmyblock);
    nmyblocks = 0;
    for(i = 0; i < varp->nblocks; i++){
        if (varp->owner[i] == nczipp->rank){
            myblocks[nmyblocks] = i;
            nmyblock++;
        }
    }

    /* 
     * Compute size to receive
     * We only need size here, packing will happen after receving
     */

    // Calculate the amount we need to receive from other process
    recvcounts = (int*)NCI_Malloc(sizeof(int) * nczipp->np);
    rdispls = (int*)NCI_Malloc(sizeof(int) * nczipp->np);
    memset(sendcounts, 0, sizeof(int) * nczipp->np);
    memset(packoff, 0, sizeof(int) * nczipp->np);
    for(i = 0; i < varp->nblocks; i++){
        if (varp->owner[i] == nczipp->rank){
            get_block_cord(varp, i, bcord);

            for(j = 0; j < nczipp->np; j++){
                // Overlapping region of this block
                get_block_overlap(varp, bcord, start_all[j], count_all[j], stride_all[j], tstart, tssize);

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
    for(i = 1; i < nczipp->np; i++){
        rdispls[i] = recvcounts[i - 1] + rdispls[i - 1];
    }

    // Allocate receive buffer
    rbuf = (char*)NCI_Malloc(rdispls[ncaipp->np - 1] + recvcounts[ncaipp->np - 1]);

    // Send the data to destination
    MPI_Alltoallv(sbuf, sendcounts, sdispls, MPI_BYTE, rbuf, recvcounts, rdispls, MPI_BYTE, nczipp->comm);

    /*
     * Next step is to pack data to block buffer
     */

    // Allocate buffer
    xbuf = (char*)NCI_Malloc(nmyblocks * bsize);

    // Main array is the whole block
    for(i = 0; i < varp->ndim; i++){
        tsize[i] = varp->stripesize[i];
    }

    // Pack data
    memset(packoff, 0, sizeof(int) * nczipp->np);
    for(i = 0; i < nmyblocks; i++){
        get_block_cord(varp, myblocks[i], bcord);

        for(j = 0; j < nczipp->np; j++){
            // Overlapping region of this block
            get_block_overlap(varp, bcord, start_all[j], count_all[j], stride_all[j], tstart, tssize);
            overlapsize = esize;
            for(k = 0; k < varp->ndim; k++){
                overlapsize *= tssize[k];
            }

            if (overlapsize > 0){
                // The block is the main array, overlapping region is the subarray
                for(k = 0; k < varp->ndim; k++){
                    tstart[k] -= bcord[k] * varp->stripesize[k];
                }
                
                // Pack type
                MPI_Type_create_subarray(varp->ndim, tsize, tssize, tstart, MPI_ORDER_C, ncmpii_nc2mpitype(xtype), &ptype);
                MPI_Type_commit(&ptype);

                // Pack data
                MPI_Unpack(rbuf + rdispls[j], overlapsize, packoff + j, xbuf + bsize * i, 1, ptype, nczipp->comm);

                // Free datatype
                MPI_Type_free(&ptype);
            }
        }
    }

    /* 
     * The buffer is now filled with data coming from all processes, it's time to compress
     */

    // compressed size and displacement
    zipsize = (int*)NCI_Malloc(sizeof(int) * nmyblock);
    zdispls = (int*)NCI_Malloc(sizeof(int) * (nmyblock + 1));
    memset(zipsize, 0, sizeof(int) * nmyblock);
    memset(zdispls, 0, sizeof(int) * (nmyblock + 1));

    // Calculate compressed data size
    for(i = 0; i < nmyblocks; i++){
        // Calculate compressed size
        // This is just estimate
        nczipp->zip->compress(xbuf + bsize * i, bsize, NULL, zipsize + i, varp->ndim, varp->stripesize, etype);
    }

    // Calculate total size
    for(i = 0; i < nmyblocks; i++){
        zdispls[0] += zipsize[i];
    }

    // Allocate buffer
    zbuf = (char*)NCI_Malloc(zdispls[0]);

    // Perform real compression
    for(i = 0; i < nmyblocks; i++){
        // Compressed the data
        // We get real size here
        nczipp->zip->compress(xbuf + bsize * i, bsize, zbuf + zdispls[i], zipsize + i, varp->ndim, varp->stripesize, etype);
        
        // Calculate offset
        zdispls[i + 1] = zdispls[i] + zipsize[i];
    }

    /*
     * Now it is time for a collective write
     * We start by syncing compressed size on all processes
     * Then, we can create variable large enough to store compressed data
     * Finally, we do collective write to store the data
     */

    // First sync on compressed block size
    // We use a all MAX reduce on all blocks
    // An alternative is to allgather and unpack the info

    // Allocate buffer
    zsize_local = (int*)NCI_Malloc(sizeof(int) * varp->nblocks);
    zsize_all = (int*)NCI_Malloc(sizeof(int) * varp->nblocks);
    zdispls_all = (int*)NCI_Malloc(sizeof(int) * varp->nblocks);
    memset(zsize_local, 0, sizeof(int) * varp->nblocks);
    memset(zsize_all, 0, sizeof(int) * varp->nblocks);
    memset(zdispls_all, 0, sizeof(int) * varp->nblocks);

    // Fill up local size
    for(i = 0; i < nmyblocks; i++){
        zsize_local[myblocks[i]] = zipsize[i];
    }

    // All reduce
    MPI_Allreduce(zsize_local, zsize_all, varp->nblocks, MPI_INT, MPI_MAX, nczipp->comm);

    // Calculate variable displacement
    zdispls_all[0] = 0;
    for(i = 1; i < varp->nblocks; i++){
        zdispls_all[i] = zsize_all[i - 1] + zdispls_all[i - 1];
    }

    // Enter redefine mode
    nczipp->driver->redef(nczipp->ncp);

    // Define dimension  for data variable
    sprintf(name, "_compressed_data_dim_%d", varp->varid)
    nczipp->driver->nczipio_def_dim(nczipp->ncp, name, zdispls_all[varp->nblocks - 1] + zsize_all[varp->nblocks - 1], &zdimid);

    // Define variable
    nczipp->driver->def_var(nczipp->ncp, name, NC_BYTE, 1, &zdimid, &(varp->datavarid));

    // Record offset in data variable
    err = nczipp->driver->put_att(nczipp->ncp, varp->varid, "_blockoffset", NC_INT, varp->nblocks, zdispls_all, MPI_INT); // Original datatype
    if (err != NC_NOERR) return err;

    // Switch to data mode
    nczipp->driver->enddef(nczipp->ncp);

    //Now, we generate a varn call to write out compressed data
    zstarts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * varp->nmyblocks);
    zcounts = (MPI_Offset**)NCI_Malloc(sizeof(MPI_Offset*) * varp->nmyblocks);
    zstarts[0] = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nmyblocks);
    zcounts[0] = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * varp->nmyblocks);
    for(i = 0; i < nmyblocks; i++){
        zstarts[i] = zstarts[0] + i;
        zcounts[i] = zcounts[0] + i;
        zstarts[i][0] = zdispls_all[myblocks[i]];
        zcounts[i][0] = zsize_all[myblocks[i]];
    }
    nczipp->driver->put_varn(nczipp->ncp, varp->datavarid, nmyblocks, zstarts, zcounts, zbuf, zdispls[nmyblocks - 1] + zsize[nmyblocks - 1], MPI_BYTE, reqMode);

    //  Free up buffers
    NCI_Free(bstart);
    NCI_Free(bend);
    NCI_Free(bcord);
    NCI_Free(tsize);
    NCI_Free(tssize);
    NCI_Free(tstart);
    NCI_Free(myblocks);
    NCI_Free(sendcounts);
    NCI_Free(sdispls);
    NCI_Free(recvcounts);
    NCI_Free(rdispls);
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
