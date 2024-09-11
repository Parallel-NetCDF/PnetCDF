/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_flush()             : dispatcher->flush()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <common.h>
#include <mpi.h>
#include <ncchkio_driver.h>
#include <pnc_debug.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include "../ncmpio/ncmpio_NC.h"
#include "ncchkio_internal.h"

int ncchkio_create (
	MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncpp) /* OUT */
{
	int err;
	int one	  = 1;
	void *ncp = NULL;
	NC_chk *ncchkp;
	PNC_driver *driver = NULL;
#ifdef PNETCDF_PROFILING
	double t0, t1;

	t0 = MPI_Wtime ();
#endif

	/* TODO: use comde to determine the true driver */
	driver = ncmpio_inq_driver ();
	if (driver == NULL) return NC_ENOTNC;

	err = driver->create (comm, path, cmode | NC_64BIT_DATA, ncid, info, &ncp);
	if (err != NC_NOERR) return err;

	/* Create a NC_chk object and save its driver pointer */
	ncchkp = (NC_chk *)NCI_Malloc (sizeof (NC_chk));
	if (ncchkp == NULL) DEBUG_RETURN_ERROR (NC_ENOMEM)

	ncchkp->path = (char *)NCI_Malloc (strlen (path) + 1);
	if (ncchkp->path == NULL) {
		NCI_Free (ncchkp);
		DEBUG_RETURN_ERROR (NC_ENOMEM)
	}
	strcpy (ncchkp->path, path);
	ncchkp->mode   = cmode | NC_WRITE;
	ncchkp->driver = driver;
	ncchkp->flag   = 0;
	ncchkp->ncp	   = ncp;
	ncchkp->comm   = comm;
	MPI_Comm_rank (comm, &(ncchkp->rank));
	MPI_Comm_size (comm, &(ncchkp->np));

	ncchkioi_init (ncchkp, 1);

	err = ncchkioi_extract_hint (ncchkp, info);
	if (err != NC_NOERR) return err;

	err = driver->put_att (ncchkp->ncp, NC_GLOBAL, "_comressed", NC_INT, 1, &one,
						   MPI_INT);  // Mark this file as compressed
	if (err != NC_NOERR) return err;

	*ncpp = ncchkp;

	// Timer array is not avaiable until init, can't use NC_CHK_TIMER_START
#ifdef PNETCDF_PROFILING
	t0 = MPI_Wtime () - t0;
	ncchkp->profile.tt[NC_CHK_TIMER_VAR_INIT] += t0;
	ncchkp->profile.tt[NC_CHK_TIMER_TOTAL] += t0;
#endif

	return NC_NOERR;
}

int ncchkio_open (
	MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncpp) {
	int err;
	int one			   = 0;
	void *ncp		   = NULL;
	NC_chk *ncchkp	   = NULL;
	PNC_driver *driver = NULL;
#ifdef PNETCDF_PROFILING
	double t0;

	t0 = MPI_Wtime ();
#endif

	/* TODO: use comde to determine the true driver */
	driver = ncmpio_inq_driver ();
	if (driver == NULL) {
		DEBUG_ASSIGN_ERROR (err, NC_ENOTNC)
		goto errout;
	}

	err = driver->open (comm, path, omode, ncid, info, &ncp);
	if (err != NC_NOERR) goto errout;

	/* Create a NC_chk object and save its driver pointer */
	ncchkp = (NC_chk *)NCI_Malloc (sizeof (NC_chk));
	if (ncchkp == NULL) {
		DEBUG_ASSIGN_ERROR (err, NC_ENOMEM)
		goto errout;
	}

	ncchkp->path = (char *)NCI_Malloc (strlen (path) + 1);
	if (ncchkp->path == NULL) {
		NCI_Free (ncchkp);
		DEBUG_ASSIGN_ERROR (err, NC_ENOMEM)
		goto errout;
	}
	strcpy (ncchkp->path, path);
	ncchkp->mode   = omode;
	ncchkp->driver = driver;
	if (ncchkp->mode & NC_WRITE) {
		ncchkp->flag = 0;
	} else {
		ncchkp->flag |= NC_MODE_RDONLY;
	}
	ncchkp->ncp	 = ncp;
	ncchkp->comm = comm;
	MPI_Comm_rank (comm, &(ncchkp->rank));
	MPI_Comm_size (comm, &(ncchkp->np));

	ncchkioi_init (ncchkp, 0);

	err = ncchkioi_extract_hint (ncchkp, info);
	if (err != NC_NOERR) goto errout;

	err = driver->get_att (ncchkp->ncp, NC_GLOBAL, "_comressed", &one,
						   MPI_INT);  // Mark this file as compressed
	if (err != NC_NOERR) {
		if (err == NC_ENOTATT) { err = NC_EINVAL; }
		goto errout;
	}

	// Not compressed file
	if (one != 1) {
		NCI_Free (ncchkp->path);
		NCI_Free (ncchkp);
		DEBUG_RETURN_ERROR (NC_EINVAL)
	}

	err = ncchkioi_get_default_chunk_dim (ncchkp);
	if (err != NC_NOERR) return err;

	ncchkioi_parse_var_info (ncchkp);

	*ncpp = ncchkp;

	// Timer array is not avaiable until init, can't use NC_CHK_TIMER_START
#ifdef PNETCDF_PROFILING
	t0 = MPI_Wtime () - t0;
	ncchkp->profile.tt[NC_CHK_TIMER_VAR_INIT] += t0;
	ncchkp->profile.tt[NC_CHK_TIMER_TOTAL] += t0;
#endif

	return NC_NOERR;

errout:
	if (ncp != NULL) { driver->close (ncchkp->ncp); }
	if (ncchkp != NULL) {
		if (ncchkp->path != NULL) { NCI_Free (ncchkp->path); }
		NCI_Free (ncchkp);
	}

	return err;
}

int ncchkio_close (void *ncdp) {
	int err = NC_NOERR;
#ifdef PNETCDF_PROFILING
	MPI_Offset put_size, get_size;
	char *_env_str = getenv ("PNETCDF_SHOW_PERFORMANCE_INFO");
#endif
	NC_chk *ncchkp = (NC_chk *)ncdp;

#ifdef PNETCDF_PROFILING
	if (_env_str != NULL && *_env_str != '0') { ncchkioi_update_statistics (ncchkp); }
#endif

	NC_CHK_TIMER_START (NC_CHK_TIMER_FINALIZE)
	NC_CHK_TIMER_START (NC_CHK_TIMER_TOTAL)

	if (ncchkp == NULL) DEBUG_RETURN_ERROR (NC_EBADID)

	if (!(ncchkp->flag & NC_MODE_RDONLY)) {
		int i;

		NC_CHK_TIMER_START (NC_CHK_TIMER_FINALIZE_META)

		err = ncchkp->driver->redef (ncchkp->ncp);
		if (err != NC_NOERR) { return err; }

		// record chunk dim
		for (i = 0; i < ncchkp->vars.cnt; i++) {
			if (ncchkp->vars.data[i].isnew) {
				err = ncchkp->driver->put_att (ncchkp->ncp, ncchkp->vars.data[i].varid, "_chunkdim",
											   NC_INT, ncchkp->vars.data[i].ndim,
											   ncchkp->vars.data[i].chunkdim, MPI_INT);
				if (err != NC_NOERR) { return err; }
				err =
					ncchkp->driver->put_att (ncchkp->ncp, ncchkp->vars.data[i].varid, "_filter",
											 NC_INT, 1, &(ncchkp->vars.data[i].filter), MPI_INT);
				if (err != NC_NOERR) { return err; }
			}
		}

		// Record recsize
		err = ncchkp->driver->put_att (ncchkp->ncp, NC_GLOBAL, "_recsize", NC_INT64, 1,
									   &(ncchkp->recsize),
									   MPI_LONG_LONG);	// Mark this file as compressed
		if (err != NC_NOERR) return err;

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_FINALIZE_META)
	}

#ifdef PNETCDF_PROFILING
	err = ncchkp->driver->inq_misc (ncchkp->ncp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
									NULL, &put_size, &get_size, NULL, NULL, NULL, NULL);
	CHK_ERR
	ncchkp->putsize += put_size;
	ncchkp->getsize += get_size;
#endif
	err = ncchkp->driver->close (ncchkp->ncp);
	CHK_ERR

	ncchkioi_cache_free (ncchkp);

	err = ncchkioi_var_list_free (&(ncchkp->vars));
	CHK_ERR

	err = ncchkioi_req_list_free (&(ncchkp->putlist));
	CHK_ERR
	err = ncchkioi_req_list_free (&(ncchkp->getlist));
	CHK_ERR

	NCI_Free (ncchkp->chunkdim);

	if (ncchkp->overlaptype != MPI_DATATYPE_NULL) { MPI_Type_free (&(ncchkp->overlaptype)); }
	if (ncchkp->max_cown_op != MPI_OP_NULL) { MPI_Op_free (&(ncchkp->max_cown_op)); }

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_FINALIZE)
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_TOTAL)

#ifdef PNETCDF_PROFILING
	if (_env_str != NULL && *_env_str != '0') {
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_PUT_SIZE,
								   (double)ncchkp->putsize / 1048576.0f);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_GET_SIZE,
								   (double)ncchkp->getsize / 1048576.0f);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_SEND_SIZE,
								   (double)ncchkp->sendsize / 1048576.0f);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_RECV_SIZE,
								   (double)ncchkp->recvsize / 1048576.0f);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NSEND, (double)ncchkp->nsend);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NRECV, (double)ncchkp->nrecv);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NREMOTE, (double)ncchkp->nremote);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NREQ, (double)ncchkp->nreq);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NLOCAL, (double)ncchkp->nlocal);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_NCHUNK, (double)ncchkp->nmychunks);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_VAR_SIZE,
								   (double)ncchkp->var_size_sum / 1048576.0f);
		ncchkioi_profile_add_time (ncchkp, NC_CHK_TIMER_VAR_ZSIZE,
								   (double)ncchkp->var_zsize_sum / 1048576.0f);

		ncchkioi_print_profile (ncchkp);
	}
#endif

err_out:;

	NCI_Free (ncchkp->path);

	NCI_Free (ncchkp);

	return err;
}

int ncchkio_enddef (void *ncdp) {
	int err = NC_NOERR, ret;
	int i;
	MPI_Offset logrecnalloc, drecnalloc;
	MPI_Offset rsize;
	NC_chk_var *varp;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	NC_CHK_TIMER_START (NC_CHK_TIMER_TOTAL)
	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT)

	drecnalloc	 = 1;
	logrecnalloc = 0;
	while (drecnalloc < ncchkp->default_recnalloc) {
		logrecnalloc++;
		drecnalloc <<= 1;
	}

	// Reserve header space
	rsize = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		varp  = ncchkp->vars.data + i;
		rsize = 0;
		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			if (varp->isrec) {
				rsize +=
					((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2;  // Atts
				rsize += ((8 + 32) + 8) * (ncchkp->default_recnalloc + logrecnalloc + 1);  // dims
				rsize += ((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) *
						 (ncchkp->default_recnalloc + 2 * logrecnalloc);  // vars
			} else {
				rsize +=
					((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2;  // Atts
				rsize += ((8 + 32) + 8) * 3;											   // dims
				rsize +=
					((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) * 3;  // vars
			}
		} else {
			rsize += ((8 + 16) + 4 + 8 + 4);  // Atts
		}
	}
	//rsize *= 2;	 // 2 times for future expension
	// Add additional reserve size
	rsize += ncchkp->hdr_reserve;

	err = ncchkp->driver->_enddef (ncchkp->ncp, rsize, 0, 0, 0);
	if (err != NC_NOERR) return err;

	err = ncchkioi_get_default_chunk_dim (ncchkp);
	if (err != NC_NOERR) return err;

	if (!(ncchkp->delay_init)) {
		int nread;
		int *lens;
		MPI_Aint *fdisps, *mdisps;
		MPI_Datatype ftype, mtype;
		MPI_Status status;
		NC_chk_var *varp;

		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

		lens   = NCI_Malloc (sizeof (int) * ncchkp->vars.cnt);
		fdisps = NCI_Malloc (sizeof (MPI_Aint) * ncchkp->vars.cnt * 2);
		mdisps = fdisps + ncchkp->vars.cnt;

		nread = 0;
		for (i = 0; i < ncchkp->vars.cnt; i++) {
			varp = ncchkp->vars.data + i;

			ncchkioi_var_init (ncchkp, varp, 0, NULL, NULL);

			if (!(varp->isnew)) {
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
			CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native",
							  MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to big endian
			ncchkioi_idx_in_swapn (varp->chunk_index, varp->nchunk + 1);
#endif
			MPI_Type_free (&ftype);
			MPI_Type_free (&mtype);
		}

		NCI_Free (lens);
		NCI_Free (fdisps);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)
	}

err_out:;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT)
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_TOTAL)

	return err;
}

int ncchkio__enddef (void *ncdp,
					 MPI_Offset h_minfree,
					 MPI_Offset v_align,
					 MPI_Offset v_minfree,
					 MPI_Offset r_align) {
	int err = NC_NOERR, ret;
	int i;
	MPI_Offset logrecnalloc, drecnalloc;
	MPI_Offset rsize;
	NC_chk_var *varp;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	NC_CHK_TIMER_START (NC_CHK_TIMER_TOTAL)
	NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT)

	drecnalloc	 = 1;
	logrecnalloc = 0;
	while (drecnalloc < ncchkp->default_recnalloc) {
		logrecnalloc++;
		drecnalloc <<= 1;
	}

	// Reserve header space
	rsize = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		varp  = ncchkp->vars.data + i;
		rsize = 0;
		if (varp->varkind == NC_CHK_VAR_COMPRESSED) {
			if (varp->isrec) {
				rsize +=
					((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2;  // Atts
				rsize += ((8 + 32) + 8) * (ncchkp->default_recnalloc + logrecnalloc + 1);  // dims
				rsize += ((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) *
						 (ncchkp->default_recnalloc + 2 * logrecnalloc);  // vars
			} else {
				rsize +=
					((8 + 16) + 4 + 8 + 4) * 8 + ((8 + 16) + 4 + 8 + 4 * varp->ndim) * 2;  // Atts
				rsize += ((8 + 32) + 8) * 3;											   // dims
				rsize +=
					((8 + 32) + 8 + 8 + (8 + 8 + (8 + 12 + 4 + 8 + 4)) + 4 + 8 + 8) * 3;  // vars
			}
		} else {
			rsize += ((8 + 16) + 4 + 8 + 4);  // Atts
		}
	}
	rsize *= 2;	 // 2 times for future expension

	err = ncchkp->driver->_enddef (ncchkp->ncp, h_minfree + rsize, v_align, v_minfree, r_align);
	if (err != NC_NOERR) return err;

	err = ncchkioi_get_default_chunk_dim (ncchkp);
	if (err != NC_NOERR) return err;

	if (!(ncchkp->delay_init)) {
		int nread;
		int *lens;
		MPI_Aint *fdisps, *mdisps;
		MPI_Datatype ftype, mtype;
		MPI_Status status;
		NC_chk_var *varp;

		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

		lens   = NCI_Malloc (sizeof (int) * ncchkp->vars.cnt);
		fdisps = NCI_Malloc (sizeof (MPI_Aint) * ncchkp->vars.cnt * 2);
		mdisps = fdisps + ncchkp->vars.cnt;

		nread = 0;
		for (i = 0; i < ncchkp->vars.cnt; i++) {
			varp = ncchkp->vars.data + i;

			err = ncchkioi_var_init (ncchkp, varp, 0, NULL, NULL);
			CHK_ERR

			if (!(varp->isnew)) {
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
			CHK_ERR_SET_VIEW (((NC *)(ncchkp->ncp))->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native",
							  MPI_INFO_NULL);

#ifdef WORDS_BIGENDIAN	// Switch back to big endian
			ncchkioi_idx_in_swapn (varp->chunk_index, varp->nchunk + 1);
#endif

			MPI_Type_free (&ftype);
			MPI_Type_free (&mtype);
		}

		NCI_Free (lens);
		NCI_Free (fdisps);

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)
	}

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT)
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_TOTAL)

err_out:;
	return err;
}

int ncchkio_redef (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->redef (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_begin_indep_data (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->begin_indep_data (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_end_indep_data (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->end_indep_data (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_abort (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	if (ncchkp == NULL) DEBUG_RETURN_ERROR (NC_EBADID)

	err = ncchkp->driver->abort (ncchkp->ncp);

	NCI_Free (ncchkp->path);
	NCI_Free (ncchkp);

	return err;
}

int ncchkio_inq (void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->inq (ncchkp->ncp, ndimsp, NULL, nattsp, xtendimp);
	if (err != NC_NOERR) return err;

	if (nvarsp != NULL) { *nvarsp = ncchkp->vars.cnt; }

	return NC_NOERR;
}

int ncchkio_inq_misc (void *ncdp,
					  int *pathlen,
					  char *path,
					  int *num_fix_varsp,
					  int *num_rec_varsp,
					  int *striping_size,
					  int *striping_count,
					  MPI_Offset *header_size,
					  MPI_Offset *header_extent,
					  MPI_Offset *recsize,
					  MPI_Offset *put_size,
					  MPI_Offset *get_size,
					  MPI_Info *info_used,
					  int *nreqs,
					  MPI_Offset *usage,
					  MPI_Offset *buf_size) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->inq_misc (ncchkp->ncp, pathlen, path, num_fix_varsp, num_rec_varsp,
									striping_size, striping_count, header_size, header_extent,
									recsize, put_size, get_size, info_used, nreqs, usage, buf_size);
	if (err != NC_NOERR) return err;

	if (num_fix_varsp != NULL) { *num_fix_varsp = ncchkp->vars.cnt; }

	if (nreqs != NULL) { *nreqs = ncchkp->putlist.nused + ncchkp->getlist.nused; }

	if (put_size != NULL) { *put_size += ncchkp->putsize; }

	if (get_size != NULL) { *get_size += ncchkp->getsize; }

	return NC_NOERR;
}

int ncchkio_cancel (void *ncdp, int num_req, int *req_ids, int *statuses) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->cancel (ncchkp->ncp, num_req, req_ids, statuses);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_wait (void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode) {
	int err = NC_NOERR, status = NC_NOERR;
	int i;
	int ncom = 0, nraw = 0;
	int *rawreqs = NULL, *comreqs = NULL;
	int *rawstats = NULL, *comstats = NULL;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	NC_CHK_TIMER_START (NC_CHK_TIMER_TOTAL)

	NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT)

	if (num_reqs < 0) {	 // NC_REQ_ALL || nreqs == NC_PUT_REQ_ALL || nreqs == NC_GET_REQ_ALL
		err = ncchkioi_wait (ncchkp, num_reqs, NULL, NULL, reqMode);
		if (status == NC_NOERR) { status = err; }
		err = ncchkp->driver->wait (ncchkp->ncp, num_reqs, NULL, NULL, reqMode);
		if (status == NC_NOERR) { status = err; }
		goto done;
	}

	if (num_reqs > 0) {
		// Count number of get and put requests
		for (i = 0; i < num_reqs; i++) {
			if (req_ids[i] & 1) { nraw++; }
		}

		// Allocate buffer
		ncom	= num_reqs - nraw;
		rawreqs = (int *)NCI_Malloc (sizeof (int) * nraw);
		CHK_PTR (rawreqs)
		comreqs = (int *)NCI_Malloc (sizeof (int) * ncom);
		CHK_PTR (comreqs)

		// Build put and get req list
		nraw = ncom = 0;
		for (i = 0; i < num_reqs; i++) {
			if (req_ids[i] & 1) {
				rawreqs[nraw++] = req_ids[i] >> 1;
			} else {
				comreqs[ncom++] = req_ids[i] >> 1;
			}
		}
	}

	if (statuses != NULL) {
		rawstats = (int *)NCI_Malloc (sizeof (int) * nraw);
		CHK_PTR (rawstats)
		comstats = (int *)NCI_Malloc (sizeof (int) * ncom);
		CHK_PTR (comstats)
	} else {
		rawstats = NULL;
		comstats = NULL;
	}

	if (nraw > 0) {
		err = ncchkp->driver->wait (ncchkp->ncp, nraw, rawreqs, rawstats, reqMode);
		if (status == NC_NOERR) { status = err; }
	}

	if (ncom > 0) {
		err = ncchkioi_wait (ncchkp, ncom, comreqs, comstats, reqMode);
		if (status == NC_NOERR) { status = err; }
	}

	// Assign stats
	if (statuses != NULL) {
		nraw = ncom = 0;
		for (i = 0; i < num_reqs; i++) {
			if (req_ids[i] & 1) {
				statuses[i] = rawstats[nraw++];
			} else {
				statuses[i] = comstats[ncom++];
			}
		}

		NCI_Free (rawstats);
		NCI_Free (comstats);
	}

	NCI_Free (rawreqs);
	NCI_Free (comreqs);

err_out:;
	if (status == NC_NOERR) status = err;
done:;

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_TOTAL)

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_WAIT)

	return status;
}

int ncchkio_set_fill (void *ncdp, int fill_mode, int *old_fill_mode) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->set_fill (ncchkp->ncp, fill_mode, old_fill_mode);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_fill_var_rec (void *ncdp, int varid, MPI_Offset recno) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->fill_var_rec (ncchkp->ncp, varid, recno);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_def_var_fill (void *ncdp, int varid, int no_fill, const void *fill_value) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->def_var_fill (ncchkp->ncp, varid, no_fill, fill_value);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_sync_numrecs (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->sync_numrecs (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_sync (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->sync (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}

int ncchkio_flush (void *ncdp) {
	int err;
	NC_chk *ncchkp = (NC_chk *)ncdp;

	err = ncchkp->driver->flush (ncchkp->ncp);
	if (err != NC_NOERR) return err;

	return NC_NOERR;
}
