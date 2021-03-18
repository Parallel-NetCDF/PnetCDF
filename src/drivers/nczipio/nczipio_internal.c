/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <common.h>
#include <mpi.h>
#include <nczipio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ncmpio/ncmpio_NC.h"
#include "nczipio_internal.h"

int nczipioi_init (NC_zip *nczipp, int isnew) {
	int err;

	nczipp->max_ndim		= 0;
	nczipp->max_chunk_size	= 0;
	nczipp->getsize			= 0;
	nczipp->putsize			= 0;
	nczipp->nmychunks		= 0;
	nczipp->nwrite			= 0;
	nczipp->cache_head		= NULL;
	nczipp->cache_tail		= NULL;
	nczipp->cache_used		= 0;
	nczipp->cache_limit		= 0;
	nczipp->cache_serial	= 0;
	nczipp->ndim			= 0;
	nczipp->chunkdim		= NULL;
	nczipp->assigned_chunks = 0;
	nczipp->cown_size		= 0;
	nczipp->max_cown_op		= MPI_OP_NULL;
	nczipp->overlaptype		= MPI_DATATYPE_NULL;

	err = nczipp->driver->inq (nczipp->ncp, NULL, NULL, NULL, &(nczipp->recdim));
	if (err != NC_NOERR) return err;

	if (isnew) {
		nczipp->recsize = 0;
	} else {
		err = nczipp->driver->get_att (nczipp->ncp, NC_GLOBAL, "_recsize", &(nczipp->recsize),
									   MPI_LONG_LONG);
		CHK_ERR	 // Mark this file as compressed
	}

	/* Initialize var list */
	err = nczipioi_var_list_init (&(nczipp->vars));
	if (err != NC_NOERR) return err;

	/* Initialize nonblocking list */
	err = nczipioi_req_list_init (&(nczipp->getlist));
	if (err != NC_NOERR) return err;
	err = nczipioi_req_list_init (&(nczipp->putlist));
	if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
	memset (&(nczipp->profile), 0, sizeof (NC_zip_timers));
	nczipp->sendsize = 0;
	nczipp->recvsize = 0;
	nczipp->nsend	 = 0;
	nczipp->nrecv	 = 0;
	nczipp->nremote	 = 0;
	nczipp->nreq	 = 0;
	nczipp->nlocal	 = 0;
#endif

err_out:;
	return err;
}

int nczipioi_parse_var_info (NC_zip *nczipp) {
	int err = NC_NOERR;
	int vid;
	int i;
	int nvar;
	int varkind;
	NC_zip_var *varp;

	int nread;
	int *lens;
	MPI_Aint *fdisps, *mdisps;
	MPI_Datatype ftype, mtype;
	MPI_Status status;

	NC_ZIP_TIMER_START (NC_ZIP_TIMER_VAR_INIT_META)

	err = nczipp->driver->inq (nczipp->ncp, NULL, &nvar, NULL, &(nczipp->recdim));
	CHK_ERR

	if (nvar > 0) {
		for (vid = 0; vid < nvar; vid++) {
			err = nczipp->driver->get_att (nczipp->ncp, vid, "_varkind", &varkind,
										   MPI_INT);  // Comressed var?
			if (err != NC_NOERR) { continue; }

			if (varkind == NC_ZIP_VAR_COMPRESSED || varkind == NC_ZIP_VAR_RAW) {
				err = nczipioi_var_list_add (&(nczipp->vars));
				if (err < 0) return err;
				varp = nczipp->vars.data + err;

				memset (varp, 0, sizeof (NC_zip_var));

				varp->varid	  = vid;
				varp->varkind = varkind;

				if (varp->varkind == NC_ZIP_VAR_COMPRESSED) {
					err = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_ndim", &(varp->ndim),
												   MPI_INT);  // Original dimensions
					if (err != NC_NOERR) return err;

					varp->dimids = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
					CHK_PTR (varp->dimids)
					varp->dimsize = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
					CHK_PTR (varp->dimsize)

					err = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_dimids",
												   varp->dimids, MPI_INT);	// Dimensiona IDs
					if (err != NC_NOERR) return err;

					for (i = 0; i < varp->ndim; i++) {
						nczipp->driver->inq_dim (nczipp->ncp, varp->dimids[i], NULL,
												 varp->dimsize + i);
					}
					if (varp->dimids[0] == nczipp->recdim) {
						varp->isrec = 1;
						if (varp->dimsize[0] < nczipp->recsize) {
							varp->dimsize[0] = nczipp->recsize;
						}
					} else {
						varp->isrec = 0;
					}

					err = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_datatype",
												   &(varp->xtype), MPI_INT);  // Original datatype
					if (err != NC_NOERR) return err;

					varp->esize	   = NC_Type_size (varp->xtype);
					varp->etype	   = ncmpii_nc2mpitype (varp->xtype);
					varp->chunkdim = NULL;
				}
			}
		}

		// Collective read index table
		if (!(nczipp->delay_init)) {
			lens = NCI_Malloc (sizeof (int) * nvar);
			CHK_PTR (lens)
			fdisps = NCI_Malloc (sizeof (MPI_Aint) * nvar * 2);
			CHK_PTR (fdisps)
			mdisps = fdisps + nvar;

			nread = 0;
			for (vid = 0; vid < nczipp->vars.cnt; vid++) {
				varp = nczipp->vars.data + vid;

				if (varp->varkind == NC_ZIP_VAR_COMPRESSED) {
					// Init var
					err = nczipioi_var_init (nczipp, varp, 0, NULL, NULL);
					CHK_ERR

					err = nczipp->driver->get_att (nczipp->ncp, varp->varid, "_metaoffset",
												   &(varp->metaoff), MPI_LONG_LONG);
					if (err == NC_NOERR) {
						lens[nread]		= sizeof (NC_zip_chunk_index_entry) * (varp->nchunk);
						fdisps[nread]	= varp->metaoff;
						mdisps[nread++] = (MPI_Aint) (varp->chunk_index);
					} else {
						varp->metaoff = -1;
						memset (varp->chunk_index, 0,
								sizeof (NC_zip_chunk_index_entry) * (varp->nchunk + 1));
					}
				}
			}

			if (nread) {
				nczipioi_sort_file_offset (nread, fdisps, mdisps, lens);

				MPI_Type_create_hindexed (nread, lens, fdisps, MPI_BYTE, &ftype);
				CHK_ERR_TYPE_COMMIT (&ftype);

				MPI_Type_create_hindexed (nread, lens, mdisps, MPI_BYTE, &mtype);
				CHK_ERR_TYPE_COMMIT (&mtype);

				// Set file view
				CHK_ERR_SET_VIEW (((NC *)(nczipp->ncp))->collective_fh,
								  ((NC *)(nczipp->ncp))->begin_var, MPI_BYTE, ftype, "native",
								  MPI_INFO_NULL);

				// Read data
				CHK_ERR_READ_AT_ALL (((NC *)(nczipp->ncp))->collective_fh, 0, MPI_BOTTOM, 1, mtype,
									 &status);

				// Restore file view
				CHK_ERR_SET_VIEW (((NC *)(nczipp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE,
								  "native", MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
				nczipioi_idx_in_swapn (varp - chunk_index, varp->nchunk + 1);
#endif

				MPI_Type_free (&ftype);
				MPI_Type_free (&mtype);
			}

			NCI_Free (lens);
			NCI_Free (fdisps);
		}
	}

	NC_ZIP_TIMER_STOP (NC_ZIP_TIMER_VAR_INIT_META)

err_out:;
	return err;
}