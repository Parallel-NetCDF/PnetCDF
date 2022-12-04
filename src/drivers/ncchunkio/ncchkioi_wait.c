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

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int ncchkioi_wait_put_reqs (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err = NC_NOERR;
	int i;
	unsigned int j;
	int nvar, nflag;
	unsigned int *flag, *flag_all;
	int *vids;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT_PUT)

	// Flag of touched vars
	nflag	 = ncchkp->vars.cnt / 32 + 1;
	flag	 = (unsigned int *)NCI_Malloc (sizeof (int) * nflag * 2);
	flag_all = flag + nflag;
	memset (flag, 0, sizeof (int) * nflag);
	for (i = 0; i < nreq; i++) {
		req = ncchkp->putlist.reqs + reqids[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}

	// Sync flag
	CHK_ERR_ALLREDUCE (flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, ncchkp->comm);

	// Build a skip list of touched vars
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) { nvar++; }
	}
	vids = (int *)NCI_Malloc (sizeof (int) * nvar);
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) { vids[nvar++] = i; }
	}

	// Perform collective buffer
	if (ncchkp->comm_unit == NC_CHK_COMM_CHUNK) {
		err = ncchkioi_iput_cb_chunk (ncchkp, nreq, reqids, stats);
	} else {
		err = ncchkioi_iput_cb_proc (ncchkp, nreq, reqids, stats);
	}
	CHK_ERR

#ifdef PNETCDF_PROFILING
	NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT_PUT_BARR)
	MPI_Barrier (ncchkp->comm);
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_WAIT_PUT_BARR)
#endif

	// Perform I/O for comrpessed variables
	err = ncchkioi_save_nvar (ncchkp, nvar, vids);
	CHK_ERR

err_out:;

	// Free buffers
	NCI_Free (vids);
	NCI_Free (flag);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_WAIT_PUT)

	return err;
}

/* Out drive currently can handle only one variable at a time
 * We pack all request as a large varn request
 */
int ncchkioi_wait_get_reqs (NC_chk *ncchkp, int nreq, int *reqids, int *stats) {
	int err = NC_NOERR;
	int i;
	unsigned int j;
	int nvar, nflag;
	unsigned int *flag, *flag_all;
	int *vids;
	NC_chk_req *req;

	NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT_GET)

	// Flag of touched vars
	nflag	 = ncchkp->vars.cnt / 32 + 1;
	flag	 = (unsigned int *)NCI_Malloc (sizeof (int) * nflag * 2);
	flag_all = flag + nflag;
	memset (flag, 0, sizeof (int) * nflag);
	for (i = 0; i < nreq; i++) {
		req = ncchkp->getlist.reqs + reqids[i];
		flag[req->varid >> 5] |= 1u << (req->varid % 32);
	}

	// Sync flag
	CHK_ERR_ALLREDUCE (flag, flag_all, nflag, MPI_UNSIGNED, MPI_BOR, ncchkp->comm);

	// Build a skip list of touched vars
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) { nvar++; }
	}
	vids = (int *)NCI_Malloc (sizeof (int) * nvar);
	nvar = 0;
	for (i = 0; i < ncchkp->vars.cnt; i++) {
		if (flag_all[i >> 5] & (1u << (i % 32))) { vids[nvar++] = i; }
	}

	// Perform I/O for comrpessed variables
	// ncchkioi_load_nvar(ncchkp, nvar, vids);

	// Perform collective buffer
	if (ncchkp->comm_unit == NC_CHK_COMM_CHUNK) {
		err = ncchkioi_iget_cb_chunk (ncchkp, nreq, reqids, stats);
	} else {
		err = ncchkioi_iget_cb_proc (ncchkp, nreq, reqids, stats);
		// ncchkioi_iget_cb_chunk(ncchkp, nreq, reqids, stats);
	}
	CHK_ERR

	NC_CHK_TIMER_START (NC_CHK_TIMER_GET_CONVERT)
	for (i = 0; i < nreq; i++) {
		req = ncchkp->getlist.reqs + reqids[i];
		if (req->buf != req->xbuf) {
			void *cbuf = (void *)req->buf;

			err = ncchkioiconvert (req->xbuf, cbuf, ncchkp->vars.data[req->varid].etype,
								   req->buftype, req->bufcount);
			CHK_ERR

			if (cbuf != req->buf) NCI_Free (cbuf);
		}
	}
	NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET_CONVERT)

err_out:;

	// Free buffers
	NCI_Free (vids);
	NCI_Free (flag);

	NC_CHK_TIMER_STOP (NC_CHK_TIMER_WAIT_GET)

	return err;
}

int ncchkioi_wait (NC_chk *ncchkp, int nreqs, int *reqids, int *stats, int reqMode) {
	int err = NC_NOERR;
	int i;
	int nput = 0, nget = 0;
	int *putreqs = NULL, *getreqs = NULL;
	int *putstats = NULL, *getstats = NULL;

	if (nreqs == NC_REQ_ALL || nreqs == NC_PUT_REQ_ALL) {
		nput	= ncchkp->putlist.nused;
		putreqs = (int *)NCI_Malloc (sizeof (int) * nput);
		CHK_PTR (putreqs)
		memcpy (putreqs, ncchkp->putlist.ids, nput * sizeof (int));
	}
	if (nreqs == NC_REQ_ALL || nreqs == NC_GET_REQ_ALL) {
		nget	= ncchkp->getlist.nused;
		getreqs = (int *)NCI_Malloc (sizeof (int) * nget);
		CHK_PTR (getreqs)
		memcpy (getreqs, ncchkp->getlist.ids, nget * sizeof (int));
	}

	if (nreqs > 0) {
		// Count number of get and put requests
		for (i = 0; i < nreqs; i++) {
			if (reqids[i] & 1) { nput++; }
		}

		// Allocate buffer
		nget	= nreqs - nput;
		putreqs = (int *)NCI_Malloc (sizeof (int) * nput);
		CHK_PTR (putreqs)
		getreqs = (int *)NCI_Malloc (sizeof (int) * nget);
		CHK_PTR (getreqs)

		// Build put and get req list
		nput = nget = 0;
		for (i = 0; i < nreqs; i++) {
			if (reqids[i] & 1) {
				putreqs[nput++] = reqids[i] >> 1;
			} else {
				getreqs[nget++] = reqids[i] >> 1;
			}
		}
	}

	if (ncchkp->delay_init) {
		NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_WAIT)
		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT)
		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_INIT_META)

		err = ncchkioi_init_nvar (ncchkp, nput, putreqs, nget, getreqs);  // nput + nget = real nreq
		if (err != NC_NOERR) { return err; }

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT)
		NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_INIT_META)
		NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT)
	} else {
		NC_CHK_TIMER_PAUSE (NC_CHK_TIMER_WAIT)
		NC_CHK_TIMER_START (NC_CHK_TIMER_VAR_RESIZE)

		// Sync number of rec
		err =
			ncchkioi_resize_nvar (ncchkp, nput, putreqs, nget, getreqs);  // nput + nget = real nreq
		if (err != NC_NOERR) { return err; }

		NC_CHK_TIMER_STOP (NC_CHK_TIMER_VAR_RESIZE)
		NC_CHK_TIMER_START (NC_CHK_TIMER_WAIT)
	}

	if (stats != NULL) {
		putstats = (int *)NCI_Malloc (sizeof (int) * nput);
		CHK_PTR (putstats)
		getstats = (int *)NCI_Malloc (sizeof (int) * nget);
		CHK_PTR (getstats)
		memset (putstats, 0, sizeof (int) * nput);
		memset (getstats, 0, sizeof (int) * nget);
	} else {
		putstats = NULL;
		getstats = NULL;
	}

	if ((ncchkp->mode & NC_WRITE) && nreqs != NC_GET_REQ_ALL) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_PUT)
		err = ncchkioi_wait_put_reqs (ncchkp, nput, putreqs, putstats);
		CHK_ERR
		NC_CHK_TIMER_STOP (NC_CHK_TIMER_PUT)
	}

	if (nreqs != NC_PUT_REQ_ALL) {
		NC_CHK_TIMER_START (NC_CHK_TIMER_GET)
		err = ncchkioi_wait_get_reqs (ncchkp, nget, getreqs, getstats);
		CHK_ERR
		NC_CHK_TIMER_STOP (NC_CHK_TIMER_GET)
	}

	// Assign stats
	if (stats != NULL) {
		nput = nget = 0;
		for (i = 0; i < nreqs; i++) {
			if (reqids[i] & 1) {
				stats[i] = putstats[nput++];
			} else {
				stats[i] = getstats[nget++];
			}
		}

		NCI_Free (putstats);
		NCI_Free (getstats);
	}

	// Remove from req list
	for (i = 0; i < nput; i++) { ncchkioi_req_list_remove (&(ncchkp->putlist), putreqs[i]); }
	for (i = 0; i < nget; i++) { ncchkioi_req_list_remove (&(ncchkp->getlist), getreqs[i]); }

err_out:;
	NCI_Free (putreqs);
	NCI_Free (getreqs);

	return err;
}
