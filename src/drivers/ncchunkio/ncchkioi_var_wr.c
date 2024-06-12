/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
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
#include <config.h>
#endif

#include <common.h>
#include <math.h>
#include <mpi.h>
#include <ncchkio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ncmpio/ncmpio_NC.h"
#include "ncchkio_internal.h"

int ncchkioi_save_var (NC_chk *ncchkp, NC_chk_var *varp) {
	int i, j, k, l, err = NC_NOERR;
	int *zsizes = NULL, *zsizes_all = NULL;
	MPI_Datatype mtype, ftype;	// Memory and file datatype
	int wcnt;
	int *lens		= NULL;
	MPI_Aint *disps = NULL;
	MPI_Status status;
	MPI_Offset *zoffs = NULL;
	MPI_Offset voff;
	void **zbufs = NULL;
	int zdimid, zvarid;
	int put_size;
	char name[128];	 // Name of objects
	NC *ncp = (NC *)(ncchkp->ncp);

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO)

	// Allocate buffer for compression
	zsizes = (int *)NCI_Malloc (sizeof (int) * varp->nchunk);
	CHK_PTR (zsizes)
	zbufs = (void **)NCI_Malloc (sizeof (void *) * varp->nmychunk);
	CHK_PTR (zbufs)
	zsizes_all = (int *)NCI_Malloc (sizeof (int) * varp->nchunk);
	CHK_PTR (zsizes_all)
	zoffs = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * (varp->nchunk + 1));
	CHK_PTR (zoffs)

	// Allocate buffer for I/O
	wcnt = 0;
	for (l = 0; l < varp->nmychunk; l++) {
		k = varp->mychunks[l];
		if (varp->dirty[k]) { wcnt++; }
	}
	if (ncchkp->rank == varp->chunk_owner[0]) { wcnt += 1; }
	lens = (int *)NCI_Malloc (sizeof (int) * wcnt);
	CHK_PTR (lens)
	disps = (MPI_Aint *)NCI_Malloc (sizeof (MPI_Aint) * wcnt);
	CHK_PTR (disps)

	memset (zsizes, 0, sizeof (int) * varp->nchunk);

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_COM)

	// Compress each chunk we own
	if (varp->filter_driver != NULL) {
		varp->filter_driver->init (MPI_INFO_NULL);
		for (l = 0; l < varp->nmychunk; l++) {
			k = varp->mychunks[l];

			if (varp->dirty[k]) {
				// Apply compression
				err = varp->filter_driver->compress_alloc (varp->chunk_cache[k]->buf, varp->chunksize,
												 zbufs + l, zsizes + k, varp->ndim, varp->chunkdim,
												 varp->etype);
				CHK_ERR
			}
		}
		varp->filter_driver->finalize ();
	} else {
		for (l = 0; l < varp->nmychunk; l++) {
			k = varp->mychunks[l];
			if (varp->dirty[k]) {
				zbufs[l]  = varp->chunk_cache[k]->buf;
				zsizes[k] = varp->chunksize;
			}
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_COM)

#ifdef PNETCDF_PROFILING
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_BARR)
	MPI_Barrier (ncchkp->comm);
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_BARR)
#endif

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_SYNC)

	// Sync compressed data size with other processes
	CHK_ERR_ALLREDUCE (zsizes, zsizes_all, varp->nchunk, MPI_INT, MPI_MAX, ncchkp->comm);

	if (varp->metaoff < 0 || varp->expanded) {
		zoffs[0] = varp->nchunkalloc * sizeof (NC_chk_chunk_index_entry);
	} else {
		zoffs[0] = 0;
	}
	for (i = 0; i < varp->nchunk; i++) { zoffs[i + 1] = zoffs[i] + zsizes_all[i]; }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_SYNC)

	if (zoffs[varp->nchunk] > 0) {	// No need to do I/O if no dirty chunk to write
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_INIT)

		/* Write comrpessed variable
		 * We start by defining data variable and writing metadata
		 * Then, we create buffer type and file type for data
		 * Finally MPI collective I/O is used for writing data
		 */

		// Enter redefine mode
		ncchkp->driver->redef (ncchkp->ncp);

		// Prepare data variable

		// Define dimension for data variable
		sprintf (name, "_datablock_dim_%d", ncchkp->nwrite);
		err = ncchkp->driver->def_dim (ncchkp->ncp, name, zoffs[varp->nchunk], &zdimid);
		if (err != NC_NOERR) return err;

		// Define data variable
		sprintf (name, "_datablock_%d", ncchkp->nwrite);
		err = ncchkp->driver->def_var (ncchkp->ncp, name, NC_BYTE, 1, &zdimid, &(zvarid));
		if (err != NC_NOERR) return err;

		// Mark as data variable
		i	= NC_CHK_VAR_DATA;
		err = ncchkp->driver->put_att (ncchkp->ncp, zvarid, "_varkind", NC_INT, 1, &i, MPI_INT);
		if (err != NC_NOERR) return err;

		// Record serial
		ncchkp->nwrite++;
		err = ncchkp->driver->put_att (ncchkp->ncp, NC_GLOBAL, "_nwrite", NC_INT, 1,
									   &(ncchkp->nwrite), MPI_INT);
		if (err != NC_NOERR) return err;

		// Metadata offset
		// Real metadata offset is only known after enddef
		// We reserve the space so we don't need to enter define mode again
		if (varp->metaoff < 0) {
			err = ncchkp->driver->put_att (ncchkp->ncp, varp->varid, "_metaoffset", NC_INT64, 1,
										   &(varp->metaoff), MPI_LONG_LONG);
			if (err != NC_NOERR) return err;
		}

		// Switch to data mode
		err = ncchkp->driver->enddef (ncchkp->ncp);
		if (err != NC_NOERR) return err;

		// Update metadata
		voff = ncp->vars.value[zvarid]->begin;
		for (i = 0; i < varp->nchunk; i++) {
			if (zsizes_all[i] > 0) {
				varp->chunk_index[i].len = zsizes_all[i];
				varp->chunk_index[i].off = zoffs[i] + voff - ncp->begin_var;
			}
		}

		if (varp->metaoff < 0 || varp->expanded) {
			varp->metaoff = voff - ncp->begin_var;
			err = ncchkp->driver->put_att (ncchkp->ncp, varp->varid, "_metaoffset", NC_INT64, 1,
										   &(varp->metaoff), MPI_LONG_LONG);
			if (err != NC_NOERR) return err;

			// unset expand flag
			varp->expanded = 0;
		}

		/* Carry out coll I/O
		 * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed
		 * when count is 0 We use a dummy call inplace of type with 0 count
		 */
		if (wcnt > 0) {
			// Create file type
			l = 0;
			if (ncchkp->rank == varp->chunk_owner[0]) {	 // First chunk owner writes metadata
				lens[l]	   = (varp->nchunk) * sizeof (NC_chk_chunk_index_entry);
				disps[l++] = (MPI_Aint)varp->metaoff + ncp->begin_var;
			}
			for (i = 0; i < varp->nmychunk; i++) {
				k = varp->mychunks[i];

				// Record compressed size
				if (varp->dirty[k]) {
					lens[l]	   = zsizes[k];
					disps[l++] = (MPI_Aint) (varp->chunk_index[k].off) + ncp->begin_var;
				}
			}
			MPI_Type_create_hindexed (wcnt, lens, disps, MPI_BYTE, &ftype);
			CHK_ERR_TYPE_COMMIT (&ftype);

			// Create memory buffer type
			l = 0;
			if (ncchkp->rank == varp->chunk_owner[0]) {	 // First chunk owner writes metadata
				lens[l]	   = (varp->nchunk) * sizeof (NC_chk_chunk_index_entry);
				disps[l++] = (MPI_Aint)varp->chunk_index;
			}
			for (i = 0; i < varp->nmychunk; i++) {
				k = varp->mychunks[i];

				// Record compressed size
				if (varp->dirty[k]) {
					lens[l]	   = zsizes[k];
					disps[l++] = (MPI_Aint)zbufs[i];
				}
			}
			err = MPI_Type_create_hindexed (wcnt, lens, disps, MPI_BYTE, &mtype);
			CHK_ERR_TYPE_COMMIT (&mtype);

			NC_CHK_TIMER_SWAP (NC_CHK_TIMER_PUT_IO_INIT, NC_CHK_TIMER_PUT_IO_WR)

#ifdef WORDS_BIGENDIAN	// NetCDF data is big endian
			if (ncchkp->rank == varp->chunk_owner[0]) {
				ncchkioi_idx_in_swapn (varp - chunk_index, varp->nchunk);
			}
#endif

			// Perform MPI-IO
			// Set file view
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
			// Write data
			CHK_ERR_WRITE_AT_ALL (ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
			// Restore file view
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
			if (ncchkp->rank == varp->chunk_owner[0]) {
				ncchkioi_idx_in_swapn (varp - chunk_index, varp->nchunk);
			}
#endif

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_WR)

#ifdef _USE_MPI_GET_COUNT
			MPI_Get_count (&status, MPI_BYTE, &put_size);
#else
			MPI_Type_size (mtype, &put_size);
#endif
			ncchkp->putsize += put_size;

			// Free type
			MPI_Type_free (&ftype);
			MPI_Type_free (&mtype);
		} else {
			NC_CHK_TIMER_SWAP (NC_CHK_TIMER_PUT_IO_INIT, NC_CHK_TIMER_PUT_IO_WR)

			// Follow coll I/O with dummy call
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
			CHK_ERR_WRITE_AT_ALL (ncp->collective_fh, 0, MPI_BOTTOM, 0, MPI_BYTE, &status);
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_WR)
		}
	}

err_out:;
	// Free buffers
	NCI_Free (zsizes);
	NCI_Free (zsizes_all);
	NCI_Free (zoffs);
	for (l = 0; l < varp->nmychunk; l++) {
		k = varp->mychunks[l];
		if (varp->dirty[k]) {
			if (varp->filter_driver != NULL) { free (zbufs[l]); }
			// Clear dirty flag
			varp->dirty[k] = 0;
		}
	}
	NCI_Free (zbufs);

	NCI_Free (lens);
	NCI_Free (disps);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO)

	return err;
}

int ncchkioi_save_nvar (NC_chk *ncchkp, int nvar, int *varids) {
	int i, j, k, l, err = NC_NOERR;
	int vid;  // Iterator for variable id
	int cid;  // Iterator for chunk id
	int total_nchunks = 0;
	int *zsizes = NULL, *zsizes_all = NULL, *zsizesp = NULL, *zsizes_allp = NULL;
	int nreq;
	MPI_Offset *zoffs = NULL, *zoffsp;
	MPI_Offset start, count, oldzoff, voff;
	MPI_Datatype mtype, ftype;	// Memory and file datatype
	int wcnt, ccnt, wcur, ccur;
	int *lens		 = NULL;
	MPI_Aint *mdisps = NULL, *fdisps = NULL;
	MPI_Status status;
	MPI_Request *reqs = NULL;
	int put_size;
	void **zbufs = NULL;
	int *zdels	 = NULL;
	int zdimid, zvarid;
	char name[128];	 // Name of objects
	NC_chk_var *varp;
	NC *ncp = (NC *)(ncchkp->ncp);
	NC_var *ncvarp;

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO)
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_INIT)

	wcnt = 0;
	ccnt = 0;
	for (i = 0; i < nvar; i++) {
		varp = ncchkp->vars.data + varids[i];
		if (ncchkp->rank == varp->chunk_owner[0]) { wcnt += 1; }
		for (l = 0; l < varp->nmychunk; l++) {
			k = varp->mychunks[l];
			if (varp->dirty[k]) { ccnt++; }
		}
		total_nchunks += varp->nchunk + 1;
	}
	wcnt += ccnt;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_INIT)

	// Allocate reqid for metadata
	reqs = (MPI_Request *)NCI_Malloc (sizeof (MPI_Request) * nvar);
	CHK_PTR (reqs)

	// Allocate buffer for compression
	zsizes = (int *)NCI_Malloc (sizeof (int) * total_nchunks);
	CHK_PTR (zsizes)
	zsizes_all = (int *)NCI_Malloc (sizeof (int) * total_nchunks);
	CHK_PTR (zsizes_all)
	zbufs = (void **)NCI_Malloc (sizeof (void *) * ccnt);
	CHK_PTR (zbufs)
	zdels = (int *)NCI_Malloc (sizeof (int) * ccnt);
	CHK_PTR (zdels)
	zoffs = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * (total_nchunks + 1));
	CHK_PTR (zoffs)

	// Allocate buffer file type
	mdisps = (MPI_Aint *)NCI_Malloc (sizeof (MPI_Aint) * wcnt);
	CHK_PTR (mdisps)
	lens = (int *)NCI_Malloc (sizeof (int) * wcnt);
	CHK_PTR (lens)
	fdisps = (MPI_Aint *)NCI_Malloc (sizeof (MPI_Aint) * wcnt);
	CHK_PTR (fdisps)

	ccur		= 0;
	zsizesp		= zsizes + nvar;
	zsizes_allp = zsizes_all + nvar;
	for (vid = 0; vid < nvar; vid++) {
		varp = ncchkp->vars.data + varids[vid];

		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_COM)

		oldzoff = zoffs[varp->nchunk];

		memset (zsizesp, 0, sizeof (int) * varp->nchunk);

		// Compress each chunk we own
		if (varp->filter_driver != NULL) {
			varp->filter_driver->init (MPI_INFO_NULL);
			for (l = 0; l < varp->nmychunk; l++) {
				cid = varp->mychunks[l];

				// Apply compression
				if (varp->dirty[cid]) {
					zdels[ccur] = 1;
					err = varp->filter_driver->compress_alloc (varp->chunk_cache[cid]->buf, varp->chunksize,
													 zbufs + (ccur++), zsizesp + cid, varp->ndim,
													 varp->chunkdim, varp->etype);
					CHK_ERR
				}
			}
			varp->filter_driver->finalize ();
		} else {
			for (l = 0; l < varp->nmychunk; l++) {
				cid = varp->mychunks[l];
				if (varp->dirty[cid]) {
					zsizesp[cid]  = varp->chunksize;
					zdels[ccur]	  = 0;
					zbufs[ccur++] = varp->chunk_cache[cid]->buf;
				}
			}
		}

		NC_CHK_TIMER_SWAP (NC_CHK_TIMER_PUT_IO_COM, NC_CHK_TIMER_PUT_IO_SYNC)

		// Sync compressed data size with other processes
		CHK_ERR_IALLREDUCE (zsizesp, zsizes_allp, varp->nchunk, MPI_INT, MPI_MAX, ncchkp->comm,
							reqs + vid);

		if (varp->metaoff < 0 || varp->expanded) {
			zsizes_all[vid] = varp->nchunkalloc * sizeof (NC_chk_chunk_index_entry);
		} else {
			zsizes_all[vid] = 0;
		}

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_SYNC)

		zsizesp += varp->nchunk;
		zsizes_allp += varp->nchunk;
	}

#ifdef PNETCDF_PROFILING
	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_BARR)
	MPI_Barrier (ncchkp->comm);
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_BARR)
#endif

	/* Write comrpessed variable
	 * We start by defining data variable and writing metadata
	 * Then, we create buffer type and file type for data
	 * Finally MPI collective I/O is used for writing data
	 */

	NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_SYNC)
	zsizes_allp = zsizes_all + nvar;
	for (vid = 0; vid < nvar; vid++) {
		varp = ncchkp->vars.data + varids[vid];
		CHK_ERR_WAIT (reqs + vid, &status);
		zsizes_allp += varp->nchunk;
	}
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_SYNC)

	zoffs[0] = 0;
	for (i = 0; i < total_nchunks; i++) { zoffs[i + 1] = zoffs[i] + zsizes_all[i]; }

	if (zoffs[total_nchunks] > 0) {	 // No need to do I/O if no dirty chunk to write
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT_IO_INIT)

		// Prepare data variable

		// Enter redefine mode
		ncchkp->driver->redef (ncchkp->ncp);

		// Define dimension for data variable
		sprintf (name, "_datablock_dim_%d", ncchkp->nwrite);
		err = ncchkp->driver->def_dim (ncchkp->ncp, name, zoffs[total_nchunks], &zdimid);
		if (err != NC_NOERR) return err;

		// Define data variable
		sprintf (name, "_datablock_%d", ncchkp->nwrite);
		err = ncchkp->driver->def_var (ncchkp->ncp, name, NC_BYTE, 1, &zdimid, &zvarid);
		if (err != NC_NOERR) return err;

		// Mark as data variable
		i	= NC_CHK_VAR_DATA;
		err = ncchkp->driver->put_att (ncchkp->ncp, zvarid, "_varkind", NC_INT, 1, &i, MPI_INT);
		if (err != NC_NOERR) return err;

		// Record serial
		ncchkp->nwrite++;
		err = ncchkp->driver->put_att (ncchkp->ncp, NC_GLOBAL, "_nwrite", NC_INT, 1,
									   &(ncchkp->nwrite), MPI_INT);
		if (err != NC_NOERR) return err;

		// Metadata offset
		for (vid = 0; vid < nvar; vid++) {
			varp = ncchkp->vars.data + varids[vid];
			// Reserve space for _metaoffset
			if (varp->metaoff < 0) {
				err = ncchkp->driver->put_att (ncchkp->ncp, varp->varid, "_metaoffset", NC_INT64, 1,
											   &(varp->metaoff), MPI_LONG_LONG);
				if (err != NC_NOERR) return err;
			}
		}

		// Switch back to data mode
		err = ncchkp->driver->enddef (ncchkp->ncp);
		if (err != NC_NOERR) return err;

		voff = ncp->vars.value[zvarid]->begin;

		wcur = ccur = 0;
		for (vid = 0; vid < nvar; vid++) {
			varp = ncchkp->vars.data + varids[vid];

			if (varp->metaoff < 0 || varp->expanded) {
				varp->metaoff = zoffs[vid] + voff - ncp->begin_var;
				err = ncchkp->driver->put_att (ncchkp->ncp, varp->varid, "_metaoffset", NC_INT64, 1,
											   &(varp->metaoff), MPI_LONG_LONG);
				if (err != NC_NOERR) return err;

				// unset expand flag
				varp->expanded = 0;
			}

			if (ncchkp->rank == varp->chunk_owner[0]) {	 // First chunk owner writes metadata
				lens[wcur]	   = varp->nchunk * sizeof (NC_chk_chunk_index_entry);
				fdisps[wcur]   = (MPI_Aint)varp->metaoff + ncp->begin_var;
				mdisps[wcur++] = (MPI_Aint) (varp->chunk_index);

				// lens[wcur] = varp->nchunk * sizeof(int);
				// fdisps[wcur] = (MPI_Aint)(varp->metaoff + ncp->begin_var + sizeof(long long) *
				// varp->nchunkalloc); mdisps[wcur++] = (MPI_Aint)(varp->data_lens);
			}
		}

		ncchkioi_sort_file_offset (wcur, fdisps, mdisps, lens);

		zsizes_allp = zsizes_all + nvar;
		zoffsp		= zoffs + nvar;
		for (vid = 0; vid < nvar; vid++) {
			varp = ncchkp->vars.data + varids[vid];

			for (cid = 0; cid < varp->nchunk; cid++) {
				if (zsizes_allp[cid] > 0) {
					varp->chunk_index[cid].len = zsizes_allp[cid];
					varp->chunk_index[cid].off = zoffsp[cid] + voff - ncp->begin_var;
				}
			}

			/* Paramemter for file and memory type
			 * We do not know variable file offset until the end of define mode
			 * We will add the displacement later
			 */
			for (i = 0; i < varp->nmychunk; i++) {
				cid = varp->mychunks[i];

				// Record parameter
				if (varp->dirty[cid]) {
					lens[wcur]	   = varp->chunk_index[cid].len;
					fdisps[wcur]   = (MPI_Aint) (varp->chunk_index[cid].off) + ncp->begin_var;
					mdisps[wcur++] = (MPI_Aint)zbufs[ccur++];
				}
			}

			// Clear dirty flag
			memset (varp->dirty, 0, varp->nchunk * sizeof (int));

			zsizes_allp += varp->nchunk;
			zoffsp += varp->nchunk;
		}

		NC_CHK_TIMER_SWAP (NC_CHK_TIMER_PUT_IO_INIT, NC_CHK_TIMER_PUT_IO_WR)

		/* Carry our coll I/O
		 * OpenMPI will fail when set view or do I/O on type created with MPI_Type_create_hindexed
		 * when count is 0 We use a dummy call inplace of type with 0 count
		 */
		if (wcnt > 0) {
			// Create file type
			MPI_Type_create_hindexed (wcnt, lens, fdisps, MPI_BYTE, &ftype);
			CHK_ERR_TYPE_COMMIT (&ftype);

			// Create memmory type
			MPI_Type_create_hindexed (wcnt, lens, mdisps, MPI_BYTE, &mtype);
			CHK_ERR_TYPE_COMMIT (&mtype);

#ifdef WORDS_BIGENDIAN	// NetCDF data is big endian
			for (vid = 0; vid < nvar; vid++) {
				varp = ncchkp->vars.data + varids[vid];
				if (ncchkp->rank == varp->chunk_owner[0]) {
					ncchkioi_idx_in_swapn (varp - chunk_index, varp->nchunk + 1);
				}
			}
#endif

			// Perform MPI-IO
			// Set file view
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
			// Write data
			CHK_ERR_WRITE_AT_ALL (ncp->collective_fh, 0, MPI_BOTTOM, 1, mtype, &status);
			// Restore file view
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
			for (vid = 0; vid < nvar; vid++) {
				varp = ncchkp->vars.data + varids[vid];
				if (ncchkp->rank == varp->chunk_owner[0]) {
					ncchkioi_idx_in_swapn (varp - chunk_index, varp->nchunk + 1);
				}
			}
#endif

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_WR)

#ifdef _USE_MPI_GET_COUNT
			MPI_Get_count (&status, MPI_BYTE, &put_size);
#else
			MPI_Type_size (mtype, &put_size);
#endif
			ncchkp->putsize += put_size;

			// Free type
			MPI_Type_free (&ftype);
			MPI_Type_free (&mtype);
		} else {
			// Follow coll I/O with dummy call
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
			CHK_ERR_WRITE_AT_ALL (ncp->collective_fh, 0, MPI_BOTTOM, 0, MPI_BYTE, &status);
			CHK_ERR_SET_VIEW (ncp->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

			NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO_WR)
		}
	}

err_out:;
	// Free buffers
	NCI_Free (zsizes);
	NCI_Free (zsizes_all);
	NCI_Free (zoffs);
	ccur = 0;
	for (i = 0; i < ccnt; i++) {
		if (zdels[i]) { free (zbufs[i]); }
	}
	NCI_Free (zbufs);
	NCI_Free (zdels);

	NCI_Free (lens);
	NCI_Free (fdisps);
	NCI_Free (mdisps);

	NCI_Free (reqs);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT_IO)

	return err;
}
