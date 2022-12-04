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
#include <ncchkio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ncmpio/ncmpio_NC.h"
#include "ncchkio_internal.h"

int ncchkioi_init (NC_chk *ncchkp, int isnew) {
	int err;

	ncchkp->max_ndim		= 0;
	ncchkp->max_chunk_size	= 0;
	ncchkp->getsize			= 0;
	ncchkp->putsize			= 0;
	ncchkp->nmychunks		= 0;
	ncchkp->nwrite			= 0;
	ncchkp->cache_head		= NULL;
	ncchkp->cache_tail		= NULL;
	ncchkp->cache_used		= 0;
	ncchkp->cache_limit		= 0;
	ncchkp->cache_serial	= 0;
	ncchkp->ndim			= 0;
	ncchkp->chunkdim		= NULL;
	ncchkp->assigned_chunks = 0;
	ncchkp->cown_size		= 0;
	ncchkp->max_cown_op		= MPI_OP_NULL;
	ncchkp->overlaptype		= MPI_DATATYPE_NULL;

	err = ncchkp->driver->inq (ncchkp->ncp, NULL, NULL, NULL, &(ncchkp->recdim));
	if (err != NC_NOERR) return err;

	if (isnew) {
		ncchkp->recsize = 0;
	} else {
		err = ncchkp->driver->get_att (ncchkp->ncp, NC_GLOBAL, "_recsize", &(ncchkp->recsize),
									   MPI_LONG_LONG);
		CHK_ERR	 // Mark this file as compressed
	}

	/* Initialize var list */
	err = ncchkioi_var_list_init (&(ncchkp->vars));
	if (err != NC_NOERR) return err;

	/* Initialize nonblocking list */
	err = ncchkioi_req_list_init (&(ncchkp->getlist));
	if (err != NC_NOERR) return err;
	err = ncchkioi_req_list_init (&(ncchkp->putlist));
	if (err != NC_NOERR) return err;

#ifdef PNETCDF_PROFILING
	memset (&(ncchkp->profile), 0, sizeof (NC_chk_timers));
	ncchkp->sendsize = 0;
	ncchkp->recvsize = 0;
	ncchkp->nsend	 = 0;
	ncchkp->nrecv	 = 0;
	ncchkp->nremote	 = 0;
	ncchkp->nreq	 = 0;
	ncchkp->nlocal	 = 0;
#endif

err_out:;
	return err;
}

int ncchkioi_parse_var_info (NC_chk *ncchkp) {
	int err = NC_NOERR, ret;
	int vid;
	int i;
	int nvar;
	int varkind;
	NC_chk_var *varp;

	int nread;
	int *lens;
	MPI_Aint *fdisps, *mdisps;
	MPI_Datatype ftype, mtype;
	MPI_Status status;

	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

	err = ncchkp->driver->inq (ncchkp->ncp, NULL, &nvar, NULL, &(ncchkp->recdim));
	CHK_ERR

	if (nvar > 0) {
		for (vid = 0; vid < nvar; vid++) {
			err = ncchkp->driver->get_att (ncchkp->ncp, vid, "_varkind", &varkind,
										   MPI_INT);  // Comressed var?
			if (err != NC_NOERR) { continue; }

			if (varkind == NC_CHK_VAR_COMPRESSED || varkind == NC_CHK_VAR_RAW) {
				err = ncchkioi_var_list_add (&(ncchkp->vars));
				if (err < 0) return err;
				varp = ncchkp->vars.data + err;

				memset (varp, 0, sizeof (NC_chk_var));

				varp->varid	  = vid;
				varp->varkind = varkind;

				if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
					err = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_ndim", &(varp->ndim),
												   MPI_INT);  // Original dimensions
					if (err != NC_NOERR) return err;

					varp->dimids = (int *)NCI_Malloc (sizeof (int) * varp->ndim);
					CHK_PTR (varp->dimids)
					varp->dimsize = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
					CHK_PTR (varp->dimsize)

					err = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_dimids",
												   varp->dimids, MPI_INT);	// Dimensiona IDs
					if (err != NC_NOERR) return err;

					for (i = 0; i < varp->ndim; i++) {
						ncchkp->driver->inq_dim (ncchkp->ncp, varp->dimids[i], NULL,
												 varp->dimsize + i);
					}
					if (varp->dimids[0] == ncchkp->recdim) {
						varp->isrec = 1;
						if (varp->dimsize[0] < ncchkp->recsize) {
							varp->dimsize[0] = ncchkp->recsize;
						}
					} else {
						varp->isrec = 0;
					}

					err = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_datatype",
												   &(varp->xtype), MPI_INT);  // Original datatype
					if (err != NC_NOERR) return err;

					varp->esize	   = NC_Type_size (varp->xtype);
					varp->etype	   = ncmpii_nc2mpitype (varp->xtype);
					varp->chunkdim = NULL;
				}
			}
		}

		// Collective read index table
		if (!(ncchkp->delay_init)) {
			lens = NCI_Malloc (sizeof (int) * nvar);
			CHK_PTR (lens)
			fdisps = NCI_Malloc (sizeof (MPI_Aint) * nvar * 2);
			CHK_PTR (fdisps)
			mdisps = fdisps + nvar;

			nread = 0;
			for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
				varp = ncchkp->vars.data + vid;

				if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
					// Init var
					err = ncchkioi_var_init (ncchkp, varp, 0, NULL, NULL);
					CHK_ERR

					ret = ncchkp->driver->get_att (ncchkp->ncp, varp->varid, "_metaoffset",
												   &(varp->metaoff), MPI_LONG_LONG);
					if (ret == NC_NOERR) {
						lens[nread]		= sizeof (NC_chk_chunk_index_entry) * (varp->nchunk);
						fdisps[nread]	= varp->metaoff;
						mdisps[nread++] = (MPI_Aint) (varp->chunk_index);
					} else {
						varp->metaoff = -1;
						memset (varp->chunk_index, 0,
								sizeof (NC_chk_chunk_index_entry) * (varp->nchunk + 1));
					}
				}
			}

			if (nread) {
				ncchkioi_sort_file_offset (nread, fdisps, mdisps, lens);

				MPI_Type_create_hindexed (nread, lens, fdisps, MPI_BYTE, &ftype);
				CHK_ERR_TYPE_COMMIT (&ftype);

				MPI_Type_create_hindexed (nread, lens, mdisps, MPI_BYTE, &mtype);
				CHK_ERR_TYPE_COMMIT (&mtype);

				// Set file view
				CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh,
								  ((NC *)(ncchkp->ncp))->begin_var, MPI_BYTE, ftype, "native",
								  MPI_INFO_NULL);

				// Read data
				CHK_ERR_READ_AT_ALL (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BOTTOM, 1, mtype,
									 &status);

				// Restore file view
				CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE,
								  "native", MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to little endian
				ncchkioi_idx_in_swapn (varp - chunk_index, varp->nchunk + 1);
#endif

				MPI_Type_free (&ftype);
				MPI_Type_free (&mtype);
			}

			for (vid = 0; vid < ncchkp->vars.cnt; vid++) {
				varp = ncchkp->vars.data + vid;
			}

			NCI_Free (lens);
			NCI_Free (fdisps);
		}
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)

err_out:;
	return err;
}