/**************************************************************************
 *
 * This file is created by Northwestern University and Argonne National
 * Laboratory
 *
 *************************************************************************/

#include "ncconfig.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#ifndef ENOERR
#define ENOERR 0
#endif
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#ifdef _MSC_VER /* Microsoft Compilers */
#include <io.h>
#else
#include <unistd.h>
#endif

#include "nc.h"
#include "ncio.h"
#include "fbits.h"
#include "rnd.h"

/* #define INSTRUMENT 1 */
#ifdef INSTRUMENT /* debugging */
#undef NDEBUG
#include <stdio.h>
#include "instr.h"
#endif

#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

#if !defined(NDEBUG) && !defined(X_INT_MAX)
#define  X_INT_MAX 2147483647
#endif

#if 0 /* !defined(NDEBUG) && !defined(X_ALIGN) */
#define  X_ALIGN 4
#else
#undef X_ALIGN
#endif

#define MAX_NC_ID 1024

static unsigned char IDalloc[MAX_NC_ID];

void
ncmpiio_free(ncio *nciop) {
  if (nciop != NULL)
    free(nciop);
}

ncio *
ncmpiio_new(const char *path, int ioflags)
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

int
ncmpiio_create(MPI_Comm comm, const char *path, int ioflags, MPI_Info info, 
            ncio **nciopp) {
  ncio *nciop;
  int i;
  int mpiomode = (MPI_MODE_RDWR | MPI_MODE_CREATE);
  int mpireturn;
  int rank;

  MPI_Comm_rank(comm, &rank);

  fSet(ioflags, NC_WRITE);

  if(path == NULL || *path == 0)
    return EINVAL;

  nciop = ncmpiio_new(path, ioflags);
  if(nciop == NULL)
    return ENOMEM;

  nciop->mpiomode = MPI_MODE_RDWR;
  nciop->mpioflags = 0;
  nciop->comm = comm;
  nciop->mpiinfo = info;

  if (fIsSet(ioflags, NC_NOCLOBBER))
    fSet(mpiomode, MPI_MODE_EXCL);
  else {
    mpireturn = MPI_File_delete((char *)path, info);
/*
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_delete error = %s\n", rank, errorString);
      return NC_EFILE;
    }
*/
  }


  mpireturn = MPI_File_open(comm, (char *)path, mpiomode, info, &nciop->collective_fh);
  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    ncmpiio_free(nciop);
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_open error = %s\n", rank, errorString);
    return NC_EOFILE;  
  }

  for (i = 0; i < MAX_NC_ID && IDalloc[i] != 0; i++);
  if (i == MAX_NC_ID) {
    ncmpiio_free(nciop);
    return NC_ENFILE;
  }
  *((int *)&nciop->fd) = i;
  IDalloc[i] = 1;

  set_NC_collectiveFh(nciop);

  *nciopp = nciop;
  return ENOERR;  
}

int
ncmpiio_open(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
          ncio **nciopp) {
  ncio *nciop;
  int i;
  int mpiomode = fIsSet(ioflags, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
  int mpireturn;
  int rank;

  MPI_Comm_rank(comm, &rank);

  if(path == NULL || *path == 0)
    return EINVAL;
 
  nciop = ncmpiio_new(path, ioflags);
  if(nciop == NULL)
    return ENOMEM;
 
  nciop->mpiomode = mpiomode;
  nciop->mpioflags = 0;
  nciop->comm = comm;
  nciop->mpiinfo = info;
 
  mpireturn = MPI_File_open(comm, (char *)path, mpiomode, info, &nciop->collective_fh);
  if (mpireturn != MPI_SUCCESS) {
    char errorString[512];
    int  errorStringLen;
    ncmpiio_free(nciop);
    MPI_Error_string(mpireturn, errorString, &errorStringLen);
    printf("%2d: MPI_File_open error = %s\n", rank, errorString);
    return NC_EOFILE;
  }
 
  for (i = 0; i < MAX_NC_ID && IDalloc[i] != 0; i++);
  if (i == MAX_NC_ID) {
    ncmpiio_free(nciop);
    return NC_ENFILE;
  }
  *((int *)&nciop->fd) = i;
  IDalloc[i] = 1;
 
  set_NC_collectiveFh(nciop);
 
  *nciopp = nciop;
  return ENOERR; 
}

int
ncmpiio_sync(ncio *nciop) {
  int mpireturn;
  int rank;

  MPI_Comm_rank(nciop->comm, &rank);

  if(NC_independentFhOpened(nciop)) {
    mpireturn = MPI_File_sync(nciop->independent_fh);
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
      return NC_EFILE;
    }
  }

 
  if(NC_collectiveFhOpened(nciop)) {
    mpireturn = MPI_File_sync(nciop->collective_fh); 
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
      return NC_EFILE;
    }
  }


  MPI_Barrier(nciop->comm);

  return ENOERR;
}

int
ncmpiio_close(ncio *nciop, int doUnlink) {
  int status = ENOERR;
  int mpireturn;
  int rank;

  MPI_Comm_rank(nciop->comm, &rank);

  if (nciop == NULL)
    return EINVAL;

  if(NC_independentFhOpened(nciop)) {
    mpireturn = MPI_File_close(&(nciop->independent_fh));
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_close error = %s\n", rank, errorString);
      return NC_EFILE;
    }
  }

 
  if(NC_collectiveFhOpened(nciop)) {
    mpireturn = MPI_File_close(&(nciop->collective_fh));  
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_close error = %s\n", rank, errorString);
      return NC_EFILE;
    }
  }

  if (doUnlink) {
    mpireturn = MPI_File_delete((char *)nciop->path, nciop->mpiinfo);
/*
    if (mpireturn != MPI_SUCCESS) {
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_delete error = %s\n", rank, errorString);
      return NC_EFILE;
    }
*/
  }

  ncmpiio_free(nciop);

  return status;
}

int
ncmpiio_move(ncio *const nciop, off_t to, off_t from, size_t nbytes) {
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
      char errorString[512];
      int  errorStringLen;
      free(buf);
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
      return NC_EREAD;
    }
    if (rank < grpsize) {
      /* read the original data @ from+movesize+rank*bufsize */
      mpireturn = MPI_File_read_at(nciop->collective_fh,
                                   from+movesize+rank*bufsize,
                                   buf, bufcount, MPI_BYTE, &mpistatus);
      if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        mpierr = 1;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_at error = %s\n", rank, errorString);
      }
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
      char errorString[512];
      int  errorStringLen;
      free(buf);
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
      return NC_EWRITE;
    }
    if (rank < grpsize) {
      /* write to new location @ to+movesize+rank*bufsize */
      mpireturn = MPI_File_write_at(nciop->collective_fh,
                                    to+movesize+rank*bufsize,
                                   buf, bufcount, MPI_BYTE, &mpistatus);
      if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        mpierr = 1;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write_at error = %s\n", rank, errorString);
      }
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
