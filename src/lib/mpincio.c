/***********************************************************************************
 *
 * This file is created by Northwestern University and Argonne National Laboratory
 *
 **********************************************************************************/

#include "ncconfig.h"
#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#ifndef ENOERR
#define ENOERR 0
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#ifdef _MSC_VER /* Microsoft Compilers */
#include <io.h>
#else
#include <unistd.h>
#endif

#include "ncio.h"
#include "fbits.h"
#include "rnd.h"

/* #define INSTRUMENT 1 */
#if INSTRUMENT /* debugging */
#undef NDEBUG
#include <stdio.h>
#include "instr.h"
#endif

#undef MIN  /* system may define MIN somewhere and complain */
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))

#if !defined(NDEBUG) && !defined(X_INT_MAX)
#define  X_INT_MAX 2147483647
#endif

#if 0 /* !defined(NDEBUG) && !defined(X_ALIGN) */
#define  X_ALIGN 4
#else
#undef X_ALIGN
#endif

#ifdef USE_MPIO /* Following code add by Jianwei Li */

#define MAX_NC_ID 1024

static unsigned char IDalloc[MAX_NC_ID];

#endif /* USE_MPIO */

#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */ 

void
ncio_free(ncio *nciop) {
  if (nciop != NULL)
    free(nciop);
}

#else /* Following interface modified as above by Jianwei Li */

static void
ncio_free(ncio *nciop)
{
	if(nciop == NULL)
		return;

	if(nciop->free != NULL)
		nciop->free(nciop->pvt);
	
	free(nciop);
}

#endif /* USE_MPIO */

#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */  

ncio *
ncio_new(const char *path, int ioflags)
{
  size_t sz_ncio = M_RNDUP(sizeof(ncio));
  size_t sz_path = M_RNDUP(strlen(path) +1); 
  ncio *nciop; 

  nciop = (ncio *) malloc(sz_ncio + sz_path);
  if (nciop == NULL) 
    return NULL;

  nciop->ioflags = ioflags;

  nciop->path = (char *) ((char *)nciop + sz_ncio);
  (void) strcpy((char *)nciop->path, path); 

  return nciop;
}

#else /* Following interface modified as above by Jianwei Li */

static ncio *
ncio_new(const char *path, int ioflags)
{
	size_t sz_ncio = M_RNDUP(sizeof(ncio));
	size_t sz_path = M_RNDUP(strlen(path) +1);
	size_t sz_ncio_pvt;
	ncio *nciop;
 
#if ALWAYS_NC_SHARE /* DEBUG */
	fSet(ioflags, NC_SHARE);
#endif

	if(fIsSet(ioflags, NC_SHARE))
		sz_ncio_pvt = sizeof(ncio_spx);
	else
		sz_ncio_pvt = sizeof(ncio_px);

	nciop = (ncio *) malloc(sz_ncio + sz_path + sz_ncio_pvt);
	if(nciop == NULL)
		return NULL;
	
	nciop->ioflags = ioflags;
	*((int *)&nciop->fd) = -1; /* cast away const */

	nciop->path = (char *) ((char *)nciop + sz_ncio);
	(void) strcpy((char *)nciop->path, path); /* cast away const */

				/* cast away const */
	*((void **)&nciop->pvt) = (void *)(nciop->path + sz_path);

	if(fIsSet(ioflags, NC_SHARE))
		ncio_spx_init(nciop);
	else
		ncio_px_init(nciop);

	return nciop;
}

#endif /* USE_MPIO */

#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
ncio_create(MPI_Comm comm, const char *path, int ioflags, MPI_Info info, 
            ncio **nciopp) {
  ncio *nciop;
  int i;
  int mpiomode = (MPI_MODE_RDWR | MPI_MODE_CREATE);
  int mpireturn, status; 


  fSet(ioflags, NC_WRITE);

  if(path == NULL || *path == 0)
    return EINVAL;

  nciop = ncio_new(path, ioflags);
  if(nciop == NULL)
    return ENOMEM;

  nciop->mpiomode = MPI_MODE_RDWR;
  nciop->mpioflags = 0;
  nciop->comm = comm;
  nciop->mpiinfo = info;

  if (fIsSet(ioflags, NC_NOCLOBBER))
    fSet(mpiomode, MPI_MODE_EXCL);
  else
    MPI_File_delete((char *)path, info);

  mpireturn = MPI_File_open(comm, (char *)path, mpiomode, info, &nciop->collective_fh);
  if (mpireturn != MPI_SUCCESS) {
    ncio_free(nciop);
    return NC_EOFILE;  
  }

  for (i = 0; i < MAX_NC_ID && IDalloc[i] != 0; i++);
  if (i == MAX_NC_ID) {
    ncio_free(nciop);
    return NC_ENFILE;
  }
  *((int *)&nciop->fd) = i;
  IDalloc[i] = 1;

  set_NC_collectiveFh(nciop);

  *nciopp = nciop;
  return ENOERR;  
}

#else /* Following interface modified as above by Jianwei Li */

int
ncio_create(const char *path, int ioflags,
	size_t initialsz,
	off_t igeto, size_t igetsz, size_t *sizehintp,
	ncio **nciopp, void **const igetvpp)
{
	ncio *nciop;
	int oflags = (O_RDWR|O_CREAT);
	int fd;
	int status;

	if(initialsz < (size_t)igeto + igetsz)
		initialsz = (size_t)igeto + igetsz;

	fSet(ioflags, NC_WRITE);

	if(path == NULL || *path == 0)
		return EINVAL;

	nciop = ncio_new(path, ioflags);
	if(nciop == NULL)
		return ENOMEM;

	if(fIsSet(ioflags, NC_NOCLOBBER))
		fSet(oflags, O_EXCL);
	else
		fSet(oflags, O_TRUNC);
#ifdef O_BINARY
	fSet(oflags, O_BINARY);
#endif
#ifdef vms
	fd = open(path, oflags, NC_DEFAULT_CREAT_MODE, "ctx=stm");
#else
	/* Should we mess with the mode based on NC_SHARE ?? */
	fd = open(path, oflags, NC_DEFAULT_CREAT_MODE);
#endif
#if 0
	(void) fprintf(stderr, "ncio_create(): path=\"%s\"\n", path);
	(void) fprintf(stderr, "ncio_create(): oflags=0x%x\n", oflags);
#endif
	if(fd < 0)
	{
		status = errno;
		goto unwind_new;
	}
	*((int *)&nciop->fd) = fd; /* cast away const */

	if(*sizehintp < NCIO_MINBLOCKSIZE || *sizehintp > NCIO_MAXBLOCKSIZE)
	{
		/* Use default */
		*sizehintp = blksize(fd);
	}
	else
	{
		*sizehintp = M_RNDUP(*sizehintp);
	}

	if(fIsSet(nciop->ioflags, NC_SHARE))
		status = ncio_spx_init2(nciop, sizehintp);
	else
		status = ncio_px_init2(nciop, sizehintp, 1);

	if(status != ENOERR)
		goto unwind_open;

	if(initialsz != 0)
	{
		status = fgrow(fd, (off_t)initialsz);
		if(status != ENOERR)
			goto unwind_open;
	}

	if(igetsz != 0)
	{
		status = nciop->get(nciop,
				igeto, igetsz,
                        	RGN_WRITE,
                        	igetvpp);
		if(status != ENOERR)
			goto unwind_open;
	}

	*nciopp = nciop;
	return ENOERR;

unwind_open:
	(void) close(fd);
	/* ?? unlink */
	/*FALLTHRU*/
unwind_new:
	ncio_free(nciop);
	return status;
}

#endif /* USE_MPIO */


#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */

int
ncio_open(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
          ncio **nciopp) {
  ncio *nciop;
  int i;
  int mpiomode = fIsSet(ioflags, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
  int mpireturn, status;
 
  if(path == NULL || *path == 0)
    return EINVAL;
 
  nciop = ncio_new(path, ioflags);
  if(nciop == NULL)
    return ENOMEM;
 
  nciop->mpiomode = mpiomode;
  nciop->mpioflags = 0;
  nciop->comm = comm;
  nciop->mpiinfo = info;
 
  mpireturn = MPI_File_open(comm, (char *)path, mpiomode, info, &nciop->collective_fh);
  if (mpireturn != MPI_SUCCESS) {
    ncio_free(nciop);
    return NC_EOFILE;
  }
 
  for (i = 0; i < MAX_NC_ID && IDalloc[i] != 0; i++);
  if (i == MAX_NC_ID) {
    ncio_free(nciop);
    return NC_ENFILE;
  }
  *((int *)&nciop->fd) = i;
  IDalloc[i] = 1;
 
  set_NC_collectiveFh(nciop);
 
  *nciopp = nciop;
  return ENOERR; 
}

#else /* Following interface modified as above by Jianwei Li */

int
ncio_open(const char *path,
	int ioflags,
	off_t igeto, size_t igetsz, size_t *sizehintp,
	ncio **nciopp, void **const igetvpp)
{
	ncio *nciop;
	int oflags = fIsSet(ioflags, NC_WRITE) ? O_RDWR : O_RDONLY;
	int fd;
	int status;

	if(path == NULL || *path == 0)
		return EINVAL;

	nciop = ncio_new(path, ioflags);
	if(nciop == NULL)
		return ENOMEM;

#ifdef O_BINARY
	fSet(oflags, O_BINARY);
#endif
#ifdef vms
	fd = open(path, oflags, 0, "ctx=stm");
#else
	fd = open(path, oflags, 0);
#endif
	if(fd < 0)
	{
		status = errno;
		goto unwind_new;
	}
	*((int *)&nciop->fd) = fd; /* cast away const */

	if(*sizehintp < NCIO_MINBLOCKSIZE || *sizehintp > NCIO_MAXBLOCKSIZE)
	{
		/* Use default */
		*sizehintp = blksize(fd);
	}
	else
	{
		*sizehintp = M_RNDUP(*sizehintp);
	}

	if(fIsSet(nciop->ioflags, NC_SHARE))
		status = ncio_spx_init2(nciop, sizehintp);
	else
		status = ncio_px_init2(nciop, sizehintp, 0);

	if(status != ENOERR)
		goto unwind_open;

	if(igetsz != 0)
	{
		status = nciop->get(nciop,
				igeto, igetsz,
                        	0,
                        	igetvpp);
		if(status != ENOERR)
			goto unwind_open;
	}

	*nciopp = nciop;
	return ENOERR;

unwind_open:
	(void) close(fd);
	/*FALLTHRU*/
unwind_new:
	ncio_free(nciop);
	return status;
}

#endif /* USE_MPIO */

#ifdef USE_MPIO /* Following interface rewritten by Jianwei Li */ 

int
ncio_sync(ncio *nciop) {
  if(NC_independentFhOpened(nciop))
    MPI_File_sync(nciop->independent_fh);
 
  if(NC_collectiveFhOpened(nciop))
    MPI_File_sync(nciop->collective_fh); 

  MPI_Barrier(nciop->comm);

  return ENOERR;
}

int
ncio_close(ncio *nciop, int doUnlink) {
  int status = ENOERR;

  if (nciop == NULL)
    return EINVAL;

  if(NC_independentFhOpened(nciop))
    MPI_File_close(&(nciop->independent_fh));
 
  if(NC_collectiveFhOpened(nciop))
    MPI_File_close(&(nciop->collective_fh));  
 
  if (doUnlink)
    MPI_File_delete((char *)nciop->path, nciop->mpiinfo);

  ncio_free(nciop);

  return status;
}

#else /* Following interface modified as above by Jianwei Li */

int 
ncio_close(ncio *nciop, int doUnlink)
{
	int status = ENOERR;

	if(nciop == NULL)
		return EINVAL;

	status = nciop->sync(nciop);

	(void) close(nciop->fd);
	
	if(doUnlink)
		(void) unlink(nciop->path);

	ncio_free(nciop);

	return status;
}

#endif /* USE_MPIO */

#ifdef USE_MPIO /* Following interface added by Jianwei Li */

int
ncio_move(ncio *const nciop, off_t to, off_t from, size_t nbytes) {
  int mpireturn, mpierr = 0, errcheck;
  const int bufsize = 4096;
  size_t movesize, bufcount;
  int rank, grpsize;
  void *buf = malloc(bufsize);
  MPI_Comm comm;
  MPI_Status mpistatus;

  comm = nciop->comm;
  MPI_Comm_size(comm, &grpsize);
  MPI_Comm_rank(comm, &rank);

  movesize = nbytes;

  while (movesize > 0) {
    /* find a proper number of processors to participate I/O */
    while (grpsize > 1 && movesize/grpsize < bufsize)
      grpsize--;
    if (grpsize > 1) {
      bufcount = bufsize;
      movesize -= bufsize*grpsize;
    } 
    else if (movesize < bufsize) {
      bufcount = movesize;
      movesize = 0;
    } 
    else {
      bufcount = bufsize;
      movesize -= bufsize;
    }

    /* reset the file view */
    mpireturn = MPI_File_set_view(nciop->collective_fh, 0, MPI_BYTE,
                                  MPI_BYTE, "native", nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS) {
      free(buf);
      return NC_EREAD;
    }
    if (rank < grpsize) {
      /* read the original data @ from+movesize+rank*bufsize */
      mpireturn = MPI_File_read_at(nciop->collective_fh,
                                   from+movesize+rank*bufsize,
                                   buf, bufcount, MPI_BYTE, &mpistatus);
      if (mpireturn != MPI_SUCCESS)
        mpierr = 1;
    }
    MPI_Allreduce(&mpierr, &errcheck, 1, MPI_INT, MPI_LOR, comm);
    if (errcheck) {
      free(buf);
      return NC_EREAD;
    }

    MPI_Barrier(comm); /* important, in case new region overlaps old region */

    /* reset the file view */
    mpireturn = MPI_File_set_view(nciop->collective_fh, 0, MPI_BYTE,
                                  MPI_BYTE, "native", nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS) {
      free(buf);
      return NC_EWRITE;
    }
    if (rank < grpsize) {
      /* write to new location @ to+movesize+rank*bufsize */
      mpireturn = MPI_File_write_at(nciop->collective_fh,
                                    to+movesize+rank*bufsize,
                                   buf, bufcount, MPI_BYTE, &mpistatus);
      if (mpireturn != MPI_SUCCESS)
        mpierr = 1;
    }
    MPI_Allreduce(&mpierr, &errcheck, 1, MPI_INT, MPI_LOR, comm);
    if (errcheck) {
      free(buf);
      return NC_EWRITE;
    }
  }

  free(buf);
  return NC_NOERR;
}

#endif /* USE_MPIO */
