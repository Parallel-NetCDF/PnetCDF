/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *	See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id$ */

#ifndef _NCIO_H_
#define _NCIO_H_

#include <stddef.h>	/* size_t */
#include <sys/types.h>	/* off_t */
#include "netcdf.h"

#include <mpi.h>

typedef struct ncio ncio;	/* forward reference */

/*
 * A value which is an invalid off_t
 */
#define OFF_NONE  ((off_t)(-1))

/*
 * Flags used by the region layer,
 *  'rflags' argument to ncio.rel() and ncio.get().
 */
#define RGN_NOLOCK	0x1	/* Don't lock region.
				 * Used when contention control handled
				 * elsewhere.
				 */
#define RGN_NOWAIT	0x2	/* return immediate if can't lock, else wait */

#define RGN_WRITE	0x4	/* we intend to modify, else read only */

#define RGN_MODIFIED	0x8	/* we did modify, else, discard */


/*
 * The next four typedefs define the signatures
 * of function pointers in struct ncio below.
 * They are not used outside of this file and ncio.h,
 * They just some casts in the ncio.c more readable.
 */
	/*
	 * Indicate that you are done with the region which begins
	 * at offset. Only reasonable flag value is RGN_MODIFIED.
	 */
typedef int ncio_relfunc(ncio *const nciop,
		 off_t offset, int rflags);

	/*
	 * Request that the region (offset, extent)
	 * be made available through *vpp.
	 */
typedef int ncio_getfunc(ncio *const nciop,
			off_t offset, size_t extent,
			int rflags,
			void **const vpp);

	/*
	 * Like memmove(), safely move possibly overlapping data.
	 * Only reasonable flag value is RGN_NOLOCK.
	 */
typedef int ncio_movefunc(ncio *const nciop, off_t to, off_t from,
			size_t nbytes, int rflags);

	/*
	 * Write out any dirty buffers to disk and
	 * ensure that next read will get data from disk.
	 */
typedef int ncio_syncfunc(ncio *const nciop);

	/*
	 * Don't call this. 
	 * Internal function called at close to
	 * free up anything hanging off pvt;
	 */
typedef void ncio_freefunc(void *const pvt);

/* Get around cplusplus "const xxx in class ncio without constructor" error */
#if defined(__cplusplus)
#define NCIO_CONST
#else
#define NCIO_CONST const
#endif

/*
 * netcdf i/o abstraction
 */
struct ncio {
	/*
	 * A copy of the ioflags argument passed in to ncio_open()
	 * or ncio_create().
	 */
	int ioflags;

	/*
	 * The file descriptor of the netcdf file.
	 * This gets handed to the user as the netcdf id.
	 */
	NCIO_CONST int fd;

	/* 
	 * The MPI File handle and the communicator 
	 * Added by Jianwei Li
	 */
#define MPI_COLLECTIVE_FH 2
#define MPI_INDEPENDENT_FH 8

#define NC_collectiveFhOpened(nciop) \
	(fIsSet(nciop->mpioflags, MPI_COLLECTIVE_FH))

#define NC_independentFhOpened(nciop) \
	(fIsSet(nciop->mpioflags, MPI_INDEPENDENT_FH))

#define set_NC_collectiveFh(nciop) \
	fSet((nciop)->mpioflags, MPI_COLLECTIVE_FH)

#define set_NC_independentFh(nciop) \
	fSet((nciop)->mpioflags, MPI_INDEPENDENT_FH)

	char* filename;
	MPI_File collective_fh;
	MPI_File independent_fh;
	int mpiomode;
	MPI_Comm comm;
	MPI_Info mpiinfo;
	int mpioflags;

	/* member functions do the work */

	ncio_relfunc *NCIO_CONST rel;

	ncio_getfunc *NCIO_CONST get;

	ncio_movefunc *NCIO_CONST move;

	ncio_syncfunc *NCIO_CONST sync;

	ncio_freefunc *NCIO_CONST free; /* Implementation private */

	/*
	 * A copy of the 'path' argument passed in to ncio_open()
	 * or ncio_create(). Used by ncabort() to remove (unlink)
	 * the file and by error messages.
	 */
	const char *path;

	/* implementation private stuff */
	void *NCIO_CONST pvt;
};

#undef NCIO_CONST

extern ncio * 
ncio_new(const char *path, int ioflags);

extern void
ncio_free(ncio *nciop);

extern int 
ncio_close(ncio *nciop, int doUnlink);

#endif /* _NCIO_H_ */
