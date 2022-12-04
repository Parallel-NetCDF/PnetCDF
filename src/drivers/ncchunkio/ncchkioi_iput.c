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

#include "ncchkio_internal.h"

static inline int ncchkioi_init_put_req (NC_chk *ncchkp,
										 NC_chk_req *req,
										 int varid,
										 const MPI_Offset *start,
										 const MPI_Offset *count,
										 const MPI_Offset *stride,
										 const void *xbuf,
										 const void *buf) {
	int err = NC_NOERR;
	int i, j, k, l;
	int *tsize, *tssize, *tstart;  // Size for sub-array type
	int overlapsize, packoff;
	MPI_Datatype ptype;	 // Pack datatype
	NC_chk_var *varp = ncchkp->vars.data + varid;

	// Zero out the request
	memset (req, 0, sizeof (NC_chk_req));

	// Record request
	req->starts	   = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *));
	req->start	   = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
	req->starts[0] = req->start;
	memcpy (req->start, start, sizeof (MPI_Offset) * varp->ndim);
	req->counts	   = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *));
	req->count	   = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
	req->counts[0] = req->count;
	memcpy (req->count, count, sizeof (MPI_Offset) * varp->ndim);
	if (stride != NULL) {
		req->stride = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim);
		memcpy (req->stride, stride, sizeof (MPI_Offset) * varp->ndim);
	}

	req->varid	  = varid;
	req->buf	  = (void *)buf;
	req->xbuf	  = (void *)xbuf;
	req->xbufs	  = (char **)NCI_Malloc (sizeof (char *));
	req->xbufs[0] = req->xbuf;
	req->nreq	  = 1;

	return err;
}

int ncchkioi_iput_var (NC_chk *ncchkp,
					   int varid,
					   const MPI_Offset *start,
					   const MPI_Offset *count,
					   const MPI_Offset *stride,
					   const void *xbuf,
					   const void *buf,
					   int *reqid) {
	int err;
	int req_id;
	NC_chk_req req;

	err = ncchkioi_init_put_req (ncchkp, &req, varid, start, count, stride, xbuf, buf);

	// Add to req list
	ncchkioi_req_list_add (&(ncchkp->putlist), &req_id);
	ncchkp->putlist.reqs[req_id] = req;

	if (reqid != NULL) { *reqid = req_id * 2 + 1; }

	return NC_NOERR;
}

static inline int ncchkioi_init_put_varn_req (NC_chk *ncchkp,
											  NC_chk_req *req,
											  int varid,
											  int nreq,
											  MPI_Offset *const *starts,
											  MPI_Offset *const *counts,
											  const void *xbuf,
											  const void *buf) {
	int err = NC_NOERR;
	int i, j;
	MPI_Offset rsize, boff;
	NC_chk_var *varp = ncchkp->vars.data + varid;

	// Zero out the request
	memset (req, 0, sizeof (NC_chk_req));

	// Record request
	req->starts = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * nreq);
	CHK_PTR (req->starts)
	req->start = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * nreq);
	CHK_PTR (req->start)
	for (i = 0; i < nreq; i++) {
		req->starts[i] = req->start + i * varp->ndim;
		memcpy (req->starts[i], starts[i], sizeof (MPI_Offset) * varp->ndim);
	}
	req->counts = (MPI_Offset **)NCI_Malloc (sizeof (MPI_Offset *) * nreq);
	CHK_PTR (req->counts)
	req->count = (MPI_Offset *)NCI_Malloc (sizeof (MPI_Offset) * varp->ndim * nreq);
	CHK_PTR (req->count)
	for (i = 0; i < nreq; i++) {
		req->counts[i] = req->count + i * varp->ndim;
		memcpy (req->counts[i], counts[i], sizeof (MPI_Offset) * varp->ndim);
	}

	// Calculate buffer for each individual request
	req->xbufs = (char **)NCI_Malloc (sizeof (char *) * nreq);
	CHK_PTR (req->xbufs)
	boff = 0;
	for (i = 0; i < nreq; i++) {
		req->xbufs[i] = (((char *)xbuf) + boff);

		// Advance pointer by size of the request
		rsize = varp->esize;
		for (j = 0; j < varp->ndim; j++) { rsize *= counts[i][j]; }
		boff += rsize;
	}

	req->varid = varid;
	req->buf   = (void *)buf;
	req->xbuf  = (void *)xbuf;
	req->nreq  = nreq;

err_out:;
	return err;
}

int ncchkioi_iput_varn (NC_chk *ncchkp,
						int varid,
						int nreq,
						MPI_Offset *const *starts,
						MPI_Offset *const *counts,
						const void *xbuf,
						const void *buf,
						int *reqid) {
	int err = NC_NOERR;
	int req_id;
	NC_chk_req req;

	if (nreq > 1) {
		err = ncchkioi_init_put_varn_req (ncchkp, &req, varid, nreq, starts, counts, xbuf, buf);
	} else {
		err = ncchkioi_init_put_req (ncchkp, &req, varid, starts[0], counts[0], NULL, xbuf, buf);
	}
	CHK_ERR

	// Add to req list
	err = ncchkioi_req_list_add (&(ncchkp->putlist), &req_id);
	CHK_ERR
	ncchkp->putlist.reqs[req_id] = req;

	if (reqid != NULL) { *reqid = req_id * 2 + 1; }

err_out:;
	return err;
}
