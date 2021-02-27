/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <common.h>
#include <math.h>
#include <nczipio_driver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nczipio_internal.h"

/* return internal size for values of specified netCDF type */
MPI_Offset NC_Type_size (nc_type type) { /* netCDF type code */
	switch (type) {
		case NC_BYTE:
			return sizeof (char);
		case NC_CHAR:
			return sizeof (char);
		case NC_SHORT:
			return sizeof (short);
		case NC_INT:
			return sizeof (int);
		case NC_FLOAT:
			return sizeof (float);
		case NC_DOUBLE:
			return sizeof (double);
		case NC_UBYTE:
			return sizeof (unsigned char);
		case NC_USHORT:
			return sizeof (unsigned short);
		case NC_UINT:
			return sizeof (unsigned int);
		case NC_INT64:
			return sizeof (long long);
		case NC_UINT64:
			return sizeof (unsigned long long);
		default:

			return 0;
	}
}

/*
 * Convert NC type to MPI type
 */
MPI_Datatype nczipioi_nc_to_mpi_type (nc_type atype) {
	switch (atype) {
		case NC_BYTE:
			return MPI_BYTE;
		case NC_CHAR:
			return MPI_CHAR;
		case NC_SHORT:
			return MPI_SHORT;
		case NC_INT:
			return MPI_INT;
		case NC_FLOAT:
			return MPI_FLOAT;
		case NC_DOUBLE:
			return MPI_DOUBLE;
	}

	return NC_NAT;
}

/*
 * Extract mpi hints and set up the flags
 */
int nczipioi_extract_hint (NC_zip *nczipp, MPI_Info info) {
	int flag;
	char value[MPI_MAX_INFO_VAL];

	// Block assignment
	MPI_Info_get (info, "nc_zip_block_mapping", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		if (strcmp (value, "static") == 0) {
			nczipp->blockmapping = NC_ZIP_MAPPING_STATIC;
		} else {
			printf ("Warning: Unknown zip method %s, using dummy\n", value);
			nczipp->blockmapping = NC_ZIP_MAPPING_STATIC;
		}
	} else {
		nczipp->blockmapping = NC_ZIP_MAPPING_STATIC;
	}

	// Messaging unit
	MPI_Info_get (info, "nc_zip_comm_unit", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		if (strcmp (value, "chunk") == 0) {
			nczipp->comm_unit = NC_ZIP_COMM_CHUNK;
		} else if (strcmp (value, "proc") == 0) {
			nczipp->comm_unit = NC_ZIP_COMM_PROC;
		} else {
			printf ("Warning: Unknown messaging unit %s, using proc\n", value);
			nczipp->comm_unit = NC_ZIP_COMM_PROC;
		}
	} else {
		nczipp->comm_unit = NC_ZIP_COMM_PROC;
	}

	// Delay init
	nczipp->delay_init = 0;
	MPI_Info_get (info, "nc_zip_delay_init", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		if (strcmp (value, "1") == 0) { nczipp->delay_init = 1; }
	}

	// Exact chunk owner assignment
	nczipp->exact_cown = 0;
	MPI_Info_get (info, "nc_zip_exact_cown", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		if (strcmp (value, "1") == 0) { nczipp->exact_cown = 1; }
	}

	// Reserve space for records
	nczipp->default_recnalloc = NC_ZIP_DEFAULT_REC_ALLOC;
	MPI_Info_get (info, "nc_zip_nrec", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) { nczipp->default_recnalloc = atoi (value); }

	// Default zipdriver
	nczipp->default_zipdriver = NC_ZIP_DRIVER_NONE;
	MPI_Info_get (info, "nc_zip_driver", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		if (strcmp (value, "none") == 0) {
			nczipp->default_zipdriver = NC_ZIP_DRIVER_NONE;
		} else if (strcmp (value, "dummy") == 0) {
			nczipp->default_zipdriver = NC_ZIP_DRIVER_DUMMY;
		} else if (strcmp (value, "zlib") == 0) {
			nczipp->default_zipdriver = NC_ZIP_DRIVER_ZLIB;
		} else if (strcmp (value, "sz") == 0) {
			nczipp->default_zipdriver = NC_ZIP_DRIVER_SZ;
		} else {
			if (nczipp->rank == 0) { printf ("Warning: Unknown zip driver %s, use none\n", value); }
		}
	}

	// Buffer size
	nczipp->cache_limit		 = 0;  // Unlimited
	nczipp->cache_limit_hint = 0;
	MPI_Info_get (info, "nc_zip_buffer_size", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) {
		sscanf (value, "%lld", &(nczipp->cache_limit_hint));

		if (nczipp->cache_limit_hint > 0) { nczipp->cache_limit = nczipp->cache_limit_hint; }
	}

	// Chunk owning size penalty
	nczipp->cown_ratio = 3;
	MPI_Info_get (info, "nc_zip_cown_ratio", MPI_MAX_INFO_VAL - 1, value, &flag);
	if (flag) { nczipp->cown_ratio = atoi (value); }

	return NC_NOERR;
}

/*
 * Export hint based on flag
 * NOTE: We only set up the hint if it is not the default setting
 *       user hint maching the default behavior will be ignored
 */
int nczipioi_export_hint (NC_zip *nczipp, MPI_Info info) {
	char value[MPI_MAX_INFO_VAL];

	MPI_Info_set (info, "nc_compression", "enable");

	switch (nczipp->blockmapping) {
		case NC_ZIP_MAPPING_STATIC:
			MPI_Info_set (info, "nc_zip_block_mapping", "static");
			break;
	}

	switch (nczipp->comm_unit) {
		case NC_ZIP_COMM_CHUNK:
			MPI_Info_set (info, "nc_zip_comm_unit", "chunk");
			break;
		case NC_ZIP_COMM_PROC:
			MPI_Info_set (info, "nc_zip_comm_unit", "proc");
			break;
	}

	// Delay inint
	if (nczipp->delay_init) {
		MPI_Info_set (info, "nc_zip_delay_init", "1");
	} else {
		MPI_Info_set (info, "nc_zip_delay_init", "0");
	}

	// Exact cown
	if (nczipp->exact_cown) {
		MPI_Info_set (info, "nc_zip_exact_cown", "1");
	} else {
		MPI_Info_set (info, "nc_zip_exact_cown", "0");
	}

	// Reserve space for records
	sprintf (value, "%lld", nczipp->default_recnalloc);
	MPI_Info_set (info, "nc_zip_nrec", value);

	// Zip driver
	switch (nczipp->default_zipdriver) {
		case NC_ZIP_DRIVER_NONE:
			MPI_Info_set (info, "nc_zip_driver", "none");
			break;
		case NC_ZIP_DRIVER_DUMMY:
			MPI_Info_set (info, "nc_zip_driver", "dummy");
			break;
		case NC_ZIP_DRIVER_ZLIB:
			MPI_Info_set (info, "nc_zip_driver", "zlib");
			break;
		case NC_ZIP_DRIVER_SZ:
			MPI_Info_set (info, "nc_zip_driver", "sz");
			break;
	}

	// Buffer size
	sprintf (value, "%lld", nczipp->cache_limit);
	MPI_Info_set (info, "nc_zip_buffer_size", value);

	return NC_NOERR;
}

int nczipioi_print_buffer_int (char *prefix, int *buf, int len) {
	int i;
	int rank, np;
	int plen, rlen;
	char *out, *outp;
	char rankstr[16];

	MPI_Comm_size (MPI_COMM_WORLD, &np);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	rlen = sprintf (rankstr, "Rank %d: ", rank);

	plen = strlen (prefix);
	out = outp = (char *)NCI_Malloc (len * 12 + 2 + plen + rlen);

	rlen = sprintf (outp, "%s ", rankstr);
	outp += rlen;
	plen = sprintf (outp, "%s ", prefix);
	outp += plen;
	for (i = 0; i < len; i++) {
		plen = sprintf (outp, "%d ", buf[i]);
		outp += plen;
	}

	printf ("%s\n", out);
	fflush (stdout);

	NCI_Free (out);

	return NC_NOERR;
}

int nczipioi_print_buffer_int64 (char *prefix, long long *buf, int len) {
	int i;
	int rank, np;
	int plen, rlen;
	char *out, *outp;
	char rankstr[16];

	MPI_Comm_size (MPI_COMM_WORLD, &np);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	rlen = sprintf (rankstr, "Rank %d: ", rank);

	plen = strlen (prefix);
	out = outp = (char *)NCI_Malloc (len * 18 + 2 + plen + rlen);

	rlen = sprintf (outp, "%s ", rankstr);
	outp += rlen;
	plen = sprintf (outp, "%s ", prefix);
	outp += plen;
	for (i = 0; i < len; i++) {
		plen = sprintf (outp, "%lld ", buf[i]);
		outp += plen;
	}

	printf ("%s\n", out);
	fflush (stdout);

	NCI_Free (out);

	return NC_NOERR;
}
#define NCZIPIOISWAP(V0, V1)  \
	fdisps[V0] ^= fdisps[V1]; \
	fdisps[V1] ^= fdisps[V0]; \
	fdisps[V0] ^= fdisps[V1]; \
	mdisps[V0] ^= mdisps[V1]; \
	mdisps[V1] ^= mdisps[V0]; \
	mdisps[V0] ^= mdisps[V1]; \
	lens[V0] ^= lens[V1];     \
	lens[V1] ^= lens[V0];     \
	lens[V0] ^= lens[V1];

void nczipioi_sort_file_offset (int len, MPI_Aint *fdisps, MPI_Aint *mdisps, int *lens) {
	int i, j, p;
	MPI_Aint at;

	if (len < 16) {
		j = 1;
		while (j) {
			j = 0;
			for (i = 0; i < len - 1; i++) {
				if (fdisps[i] > fdisps[i + 1]) {
					NCZIPIOISWAP (i, i + 1);
					j = 1;
				}
			}
		}
	} else {
		j = len / 2;
		p = len - 1;
		NCZIPIOISWAP (j, p);

		for (i = j = 0; i < len; i++) {
			if (fdisps[i] < fdisps[p]) {
				if (i != j) { NCZIPIOISWAP (i, j); }
				j++;
			}
		}

		NCZIPIOISWAP (p, j);

		nczipioi_sort_file_offset (j, fdisps, mdisps, lens);
		nczipioi_sort_file_offset (len - j - 1, fdisps + j + 1, mdisps + j + 1, lens + j + 1);
	}
}

int nczipioi_subarray_off_len (int ndim, int *tsize, int *tssize, int *tstart, int *off, int *len) {
	int err;
	int i;

	// Try single row
	err = 0;
	for (i = 0; i < ndim - 1; i++) {
		if (tssize[i] > 1) {
			err = -1;
			break;
		}
	}
	if (err) {
		// Try contiguous block
		err = 0;
		for (i = 1; i < ndim; i++) {
			if (tssize[i] < tsize[i]) {
				err = -1;
				break;
			}
		}
		if (!err) {
			*len = 1;
			for (i = 0; i < ndim; i++) { (*len) *= tssize[i]; }
		}
	} else {
		*len = tssize[ndim - 1];
	}

	if (!err) {
		*off = 0;
		for (i = 0; i < ndim; i++) { (*off) = (*off) * tsize[i] + tstart[i]; }
	}

	return err;
}

#ifdef PNETCDF_PROFILING
int nczipioi_update_statistics (NC_zip *nczipp) {
	int i, j;
	int cid;
	NC_zip_var *varp;

	nczipp->var_size_sum = nczipp->var_zsize_sum = 0;
	for (i = 0; i < nczipp->vars.cnt; i++) {
		varp = nczipp->vars.data + i;
		if (varp->varkind == NC_ZIP_VAR_COMPRESSED) {
			for (j = 0; j < varp->nmychunk; j++) {
				cid = varp->mychunks[j];
				nczipp->var_zsize_sum += varp->chunk_index[cid].len;
			}
			nczipp->var_size_sum += varp->nmychunk * varp->chunksize;
		}
	}

	return NC_NOERR;
}
#endif

int nczipioi_get_default_chunk_dim (NC_zip *nczipp) {
	int err;
	int i;
	int ndim, dimid;
	int len;
	char *cur, *pre;
	char name[1024];
	char *env = getenv ("PNETCDF_DEFAULT_CHUNK_DIM");

	if (env != NULL) {
		err = nczipp->driver->inq (nczipp->ncp, &ndim, NULL, NULL, NULL);
		if (err != NC_NOERR) return err;

		if (ndim > nczipp->ndim) {
			nczipp->chunkdim = NCI_Realloc (nczipp->chunkdim, ndim * sizeof (int));
			for (i = nczipp->ndim; i < ndim; i++) { nczipp->chunkdim[i] = 0; }
			nczipp->ndim = ndim;
		}

		cur = pre = env;
		for (cur = pre = env; (*cur) != '\0'; cur++) {
			if ((*cur) == ';') {
				if (sscanf (pre, "%s : %d ;", name, &len) == 2) {
					if (len > 0) {
						err = nczipp->driver->inq_dimid (nczipp->ncp, name, &dimid);
						if (err == NC_NOERR) { nczipp->chunkdim[dimid] = len; }
					}
				}
				pre = cur + 1;
			}
		}
	}

	return NC_NOERR;
}

/* in-place byte swap */
void nczipioi_idx_in_swapn (NC_zip_chunk_index_entry *idx, MPI_Offset nelems) {
	NC_zip_chunk_index_entry *bufp;

	for (bufp = idx; bufp < idx + nelems; bufp++) {
		bufp->off = ((bufp->off & 0x00000000000000FFULL) << 56) |
					((bufp->off & 0x000000000000FF00ULL) << 40) |
					((bufp->off & 0x0000000000FF0000ULL) << 24) |
					((bufp->off & 0x00000000FF000000ULL) << 8) |
					((bufp->off & 0x000000FF00000000ULL) >> 8) |
					((bufp->off & 0x0000FF0000000000ULL) >> 24) |
					((bufp->off & 0x00FF000000000000ULL) >> 40) |
					((bufp->off & 0xFF00000000000000ULL) >> 56);
		bufp->len = ((bufp->len) << 24) | (((bufp->len) & 0x0000ff00) << 8) |
					(((bufp->len) & 0x00ff0000) >> 8) | (((bufp->len) >> 24));
	}
}