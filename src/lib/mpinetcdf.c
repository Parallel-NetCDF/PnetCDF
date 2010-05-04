/*********************************************************************************
 * 
 * This file is created by Northwestern University and Argonne National 
 * Laboratory
 *
 ********************************************************************************/

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>
#include <arpa/inet.h>   /* htonl(), htons() */

#include "ncmpidtype.h"

/* Local prototypes */
static int length_of_mpitype(MPI_Datatype);

const char *
ncmpi_inq_libvers(void) {
	return "version = " PNETCDF_VERSION "of 02 November 2009";
}

/* Prototypes for functions used only in this file */
static int echar(nc_type nctype,MPI_Datatype mpitype);
static int need_convert(nc_type nctype,MPI_Datatype mpitype);
static int need_swap(nc_type nctype,MPI_Datatype mpitype);
static int x_putn_schar(void *xbuf, const void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_putn_short(void *xbuf, const void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_putn_int(void *xbuf, const void *buf, MPI_Offset nelems,
		      MPI_Datatype datatype);
static int x_putn_float(void *xbuf, const void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_putn_double(void *xbuf, const void *buf, MPI_Offset nelems,
			 MPI_Datatype datatype);
static int x_getn_schar(const void *xbuf, void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_getn_short(const void *xbuf, void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_getn_int(const void *xbuf, void *buf, MPI_Offset nelems,
		      MPI_Datatype datatype);
static int x_getn_float(const void *xbuf, void *buf, MPI_Offset nelems,
			MPI_Datatype datatype);
static int x_getn_double(const void *xbuf, void *buf, MPI_Offset nelems,
			 MPI_Datatype datatype);
static int set_var_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp);
static int set_vara_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp,
			     const MPI_Offset start[], const MPI_Offset count[],
			     int getnotput);
static int set_vars_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, 
			     const MPI_Offset start[], const MPI_Offset count[], 
			     const MPI_Offset stride[], int getnotput);
static int check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm,
		       int collective);
static int check_recsize_too_big(NC *ncp);

static int ncmpi_coll_wait(NCMPI_Request request);
static int ncmpi_coll_waitall(int count, NCMPI_Request array_of_requests[]);

/* Begin Of Dataset Functions */

int 
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    int status = NC_NOERR;
    MPI_Offset sizeof_off_t = 0;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;	/* might be a good thing to hint later */
    NC *ncp;

    ncp = ncmpii_new_NC(&chunksize); /* allocate buffer for header */
    if (ncp == NULL) 
        return NC_ENOMEM;

    assert(ncp->flags == 0);

    if (fIsSet(cmode, NC_64BIT_OFFSET)) {
        /* unlike serial netcdf, we will not bother to support
         * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
         * serial netcdf has proven it's possible if datasets are small, but
         * that's a hassle we don't want to worry about */
        if (sizeof(off_t) != 8)
            return NC_ESMALL;
        fSet(ncp->flags, NC_64BIT_OFFSET);
        sizeof_off_t = 8;
    } else if (fIsSet(cmode, NC_64BIT_DATA)) {
        if (sizeof(MPI_Offset) <  8)
            return NC_ESMALL;
        fSet(ncp->flags, NC_64BIT_DATA);
        sizeof_off_t = 8;
    } else {
        fSet(ncp->flags, NC_32BIT);
        sizeof_off_t = 4;
    }
    assert(ncp->xsz = ncmpii_hdr_len_NC(ncp, sizeof_off_t));

    fSet(ncp->flags, NC_NOFILL);

    status = ncmpiio_create(comm, path, cmode, info, &ncp->nciop);  
    if (status != NC_NOERR) {
        ncmpii_free_NC(ncp);
        return status;
    }

    fSet(ncp->flags, NC_CREAT);

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
        /*
         * NC_SHARE implies sync up the number of records as well.
         * (File format version one.)
         * Note that other header changes are not shared
         * automatically.  Some sort of IPC (external to this package)
         * would be used to trigger a call to ncmpi_sync().
         */ 
        fSet(ncp->flags, NC_NSYNC);  /* sync numrecs */
        fSet(ncp->flags, NC_HSYNC);  /* sync header */
    }

    ncmpii_add_to_NCList(ncp);
    *ncidp = ncp->nciop->fd;

    return status;
}

int
ncmpi_open(MPI_Comm    comm,
           const char *path,
           int         omode,
           MPI_Info    info,
           int        *ncidp)
{
    int status = NC_NOERR;
    NC *ncp;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;	/* might be a good thing to hint later */
  
    ncp = ncmpii_new_NC(&chunksize);
    if (ncp == NULL)
        return NC_ENOMEM;

    status = ncmpiio_open(comm, path, omode, info, &ncp->nciop);
    if (status != NC_NOERR) {
        ncmpii_free_NC(ncp);
        return status;
    } 

    assert(ncp->flags == 0); 

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
        /*
         * NC_SHARE implies sync up the number of records as well.
         * (File format version one.)
         * Note that other header changes are not shared
         * automatically.  Some sort of IPC (external to this package)
         * would be used to trigger a call to ncmpi_sync().
         */ 
        fSet(ncp->flags, NC_NSYNC);  /* sync numrecs */
        fSet(ncp->flags, NC_HSYNC);  /* sync header */
    }

    status = ncmpii_hdr_get_NC(ncp); /* read header from file */
    if (status != NC_NOERR) {
        ncmpiio_close(ncp->nciop, 0);
        ncmpii_free_NC(ncp);
        return status;
    }

    ncmpii_add_to_NCList(ncp);
    *ncidp = ncp->nciop->fd;

    return status;
}

int
ncmpi_inq_format(int ncid, int *formatp)
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if(status != NC_NOERR)
        return status;

    if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
        *formatp = NC_FORMAT_64BIT_DATA;
    } else if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) {
        *formatp = NC_FORMAT_64BIT;
    } else if (fIsSet(ncp->flags, NC_32BIT)){
          *formatp = NC_FORMAT_CLASSIC;
    } else {
          *formatp = NC_FORMAT_UNKNOWN;
    }
    return 0;
}

int
ncmpi_inq_file_format(char *filename, int *formatp)
{
    int status;
    NC *ncp;
    int ncid;
	

    status = ncmpi_open(MPI_COMM_SELF, filename, 0, MPI_INFO_NULL, &ncid);
    if (status == NC_ENOTNC){ 
        *formatp = NC_FORMAT_UNKNOWN;
    }
    if(status != NC_NOERR)
        return status;
    status = ncmpii_NC_check_id(ncid, &ncp);
    if(status != NC_NOERR)
         return status;


    if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
        *formatp = NC_FORMAT_64BIT_DATA;
    } else if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) {
        *formatp = NC_FORMAT_64BIT;
    } else if (fIsSet(ncp->flags, NC_32BIT)){
        *formatp = NC_FORMAT_CLASSIC;
    } else {
        *formatp = NC_FORMAT_UNKNOWN;
    }
    status = ncmpi_close(ncid);
       
    return 0;
}


int
ncmpi_get_file_info(int ncid, MPI_Info *info_used) {
  int status = NC_NOERR;
  int mpireturn;
  NC *ncp;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if (status != NC_NOERR)
    return status;

  mpireturn = MPI_File_get_info(ncp->nciop->collective_fh, info_used);
  if (mpireturn != MPI_SUCCESS) {
      int rank;
      MPI_Comm_rank(ncp->nciop->comm, &rank);
      ncmpii_handle_error(rank, mpireturn, "MPI_File_get_info");
      return NC_EFILE;
  }
  return status;
}

int
ncmpi_redef(int ncid) {
    int status;
    NC *ncp;
    MPI_Offset mynumrecs, numrecs;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) 
        return status; 

    if (NC_readonly(ncp)) 
        return NC_EPERM;

    if (NC_indef(ncp))
        return NC_EINDEFINE;
 
    /* ensure exiting define mode always entering collective data mode */
    if (NC_indep(ncp))
        ncmpi_end_indep_data(ncid);

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
        /* re-read the header from file */
        status = ncmpii_read_NC(ncp);
        if (status != NC_NOERR)
            return status;
    } else {
        /* before enter define mode, the number of records may increase by
           independent APIs, i.e. ncp->numrecs may be incoherent and need
           to sync across all processes 
           Note that only ncp->numrecs in the header can be incoherent.
         */
        mynumrecs = ncp->numrecs;
        MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX, ncp->nciop->comm);
        if (numrecs > ncp->numrecs) {
            ncp->numrecs = numrecs;
            set_NC_ndirty(ncp);
        }
    }

    ncp->old = ncmpii_dup_NC(ncp);
    if (ncp->old == NULL)
        return NC_ENOMEM;

    fSet(ncp->flags, NC_INDEF);

    return NC_NOERR;
}

static
int update_numrecs(NC         *ncp,
                   MPI_Offset  newnumrecs)
{
    int status = NC_NOERR;
    /* update the number of records in NC and write to file header, if
       necessary */
    if (ncp->numrecs < newnumrecs) {
        ncp->numrecs = newnumrecs;
        set_NC_ndirty(ncp);
    }
    if (NC_doNsync(ncp)) {
        int      localChange=0, doChange;
        MPI_Comm comm = ncp->nciop->comm;

        if (NC_ndirty(ncp)) localChange = 1;
        MPI_Allreduce(&localChange, &doChange, 1, MPI_INT, MPI_MAX, comm);

        if (doChange) {
            /* all proc must agree on numrecs because this func is collective */
            MPI_Allreduce(&newnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX, comm);
            status = ncmpii_write_numrecs(ncp);
            if (status != NC_NOERR)
                return status;

            /* fsync to disk */
            ncmpiio_sync(ncp->nciop);
        }
    }
    return status;
}
 
int
ncmpi_begin_indep_data(int ncid) {
    int mpireturn, status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (NC_indef(ncp))  /* must not be in define mode */
        return NC_EINDEFINE;

    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_EINDEP; /* Should we skip this error? */
 
    if (!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
        /* do memory and file sync for numrecs */
        MPI_Offset newnumrecs = ncp->numrecs;
        MPI_Allreduce(&newnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX, ncp->nciop->comm);
        status = ncmpii_write_numrecs(ncp);
        if (status != NC_NOERR)
            return status;

        mpireturn = MPI_File_sync(ncp->nciop->collective_fh);   /* collective */
        if (mpireturn != MPI_SUCCESS) {
            int rank;
            MPI_Comm_rank(ncp->nciop->comm, &rank);
	    ncmpii_handle_error(rank, mpireturn, "MPI_File_sync");
            MPI_Finalize();
            return NC_EFILE;
        }
    }

    fSet(ncp->flags, NC_INDEP);

    MPI_Barrier(ncp->nciop->comm);

    return status;  
}

int 
ncmpi_end_indep_data(int ncid) {
    int mpireturn, status = NC_NOERR;
    NC *ncp;
 
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status; 

    if (!NC_indep(ncp))
        return NC_ENOTINDEP;

    if (!NC_readonly(ncp)) {
        MPI_Offset newnumrecs = ncp->numrecs;
        MPI_Allreduce(&newnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX, ncp->nciop->comm);
        status = ncmpii_write_numrecs(ncp);
        if (status != NC_NOERR)
            return status;

        /* calling file sync for those already open the file */
        if (NC_independentFhOpened(ncp->nciop)) {
            mpireturn = MPI_File_sync(ncp->nciop->independent_fh); /* independent */
            if (mpireturn != MPI_SUCCESS) {
                int rank;
                MPI_Comm_rank(ncp->nciop->comm, &rank);
	        ncmpii_handle_error(rank, mpireturn, "MPI_File_sync");
                MPI_Finalize();
                return NC_EFILE;
            }
        }
    }

    fClr(ncp->flags, NC_INDEP);
 
    MPI_Barrier(ncp->nciop->comm);

    return status;
}

int
ncmpi_enddef(int ncid) {
    int status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR)
        return status;

    if (!NC_indef(ncp)) /* must currently in define mode */
        return(NC_ENOTINDEFINE);

    return ncmpii_NC_enddef(ncp);
}

int
ncmpi_sync(int ncid) {
    int status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (NC_indef(ncp)) 
        return NC_EINDEFINE;

    if (NC_readonly(ncp))
        return ncmpii_read_NC(ncp);

    status = ncmpii_NC_sync(ncp, 0); /* write header to file */
    if (status != NC_NOERR)
        return status;

    /* calling MPI_File_sync() */
    return ncmpiio_sync(ncp->nciop);
}

int
ncmpi_abort(int ncid) {
   /*
    * In data mode, same as ncmpiio_close.
    * In define mode, descard new definition.
    * In create, remove the file.
    */
    int status, doUnlink = 0;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if(status != NC_NOERR)
        return status;

    doUnlink = NC_IsNew(ncp);

    if (ncp->old != NULL) {
        /* a plain redef, not a create */
        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_INDEF));
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
        fClr(ncp->flags, NC_INDEF);
    } 
    else if (!NC_readonly(ncp) && !NC_indef(ncp)) {
        /* data mode, write */
        status = ncmpii_NC_sync(ncp, 0);
        if (status != NC_NOERR)
            return status;
    }

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
        /* calling MPI_File_sync() */
        ncmpiio_sync(ncp->nciop);
    }
    ncmpiio_close(ncp->nciop, doUnlink);
    ncp->nciop = NULL;

    ncmpii_del_from_NCList(ncp);

    ncmpii_free_NC(ncp);

    return NC_NOERR;
}

int
ncmpi_close(int ncid) {
    int status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    /* release NC object, close the file and write dirty numrecs if necessary */
    return ncmpii_NC_close(ncp);
}

/* ncmpi_delete:
 * doesn't do anything to release resources, so call ncmpi_close before calling
 * this function.
 *
 * filename: the name of the
 * file we will remove.  info: mpi info, in case underlying file system needs
 * hints.
 */
int
ncmpi_delete(char *filename, MPI_Info info)
{
	int status = NC_NOERR;
	status = MPI_File_delete(filename, MPI_INFO_NULL);
	if (status != MPI_SUCCESS)
		return NC_EFILE;
	return NC_NOERR;
}

/* ncmpi_set_fill:
 * not actually implemented.  Anything other than NC_NOFILL is not supported.
 * Many codes use NC_NOFILL anyway, so this just gets us more source-portable
 * with existings serial netcdf codes.   Also provides a placeholder if someday
 * someone wants to implement all of set_fill 
 */
int
ncmpi_set_fill(int ncid, int fillmode, int *old_mode_ptr)
{
	int status = NC_NOERR;
	if (fillmode != NC_NOFILL)
		status = NC_EINVAL;
	return status;

}
		

/* End Of Dataset Functions */

/*
 *  MAPPING:  MPI DATATYPE   <--->   NETCDF DATATYPE	   DESCRIPTION
 *		MPI_UNSIGNED_CHAR	NC_BYTE		  uchar integer
 *		MPI_BYTE		NC_BYTE		  schar integer
 *		MPI_CHAR		NC_CHAR		  char(text)
 *		MPI_SHORT		NC_SHORT	  short
 *		MPI_INT			NC_INT		  int
 *		MPI_FLOAT		NC_FLOAT	  float
 *		MPI_DOUBLE		NC_DOUBLE	  double
 *
 *
 *  Assume: MPI_Datatype and nc_type are both enumerable types
 */

static int
length_of_mpitype(MPI_Datatype datatype) {
  if ( datatype == MPI_UNSIGNED_CHAR || 
       datatype == MPI_BYTE || 
       datatype == MPI_CHAR)
    return (int) sizeof(char);
  else if (datatype == MPI_SHORT)  return (int) sizeof(short);
  else if (datatype == MPI_INT)    return (int) sizeof(int);
  else if (datatype == MPI_LONG)   return (int) sizeof(long);
  else if (datatype == MPI_FLOAT)  return (int) sizeof(float);
  else if (datatype == MPI_DOUBLE) return (int) sizeof(double);
  else 
	fprintf(stderr, "FIXME: unknown type passed to length_of_mpitype\n");

  return -1;
}

static int
echar(nc_type nctype,MPI_Datatype mpitype) {
  return ((nctype == NC_CHAR) == (mpitype != MPI_CHAR));
}

static int
need_convert(nc_type nctype,MPI_Datatype mpitype) {
  return !( (nctype == NC_CHAR && mpitype == MPI_CHAR) ||
            (nctype == NC_BYTE && mpitype == MPI_BYTE) ||
	    (nctype == NC_BYTE && mpitype == MPI_UNSIGNED_CHAR) ||
            (nctype == NC_SHORT && mpitype == MPI_SHORT) ||
            (nctype == NC_INT && mpitype == MPI_INT) ||
	    (nctype == NC_INT && mpitype == MPI_LONG && X_SIZEOF_INT == SIZEOF_LONG) ||
            (nctype == NC_FLOAT && mpitype == MPI_FLOAT) ||
            (nctype == NC_DOUBLE && mpitype == MPI_DOUBLE) );
}

static int 
need_swap(nc_type nctype,MPI_Datatype mpitype) {
#ifdef WORDS_BIGENDIAN
  return 0;
#else
  return ( (nctype == NC_SHORT && mpitype == MPI_SHORT) ||
          (nctype == NC_INT && mpitype == MPI_INT) ||
	  (nctype == NC_INT && mpitype == MPI_LONG && X_SIZEOF_INT == SIZEOF_LONG) ||
          (nctype == NC_FLOAT && mpitype == MPI_FLOAT) ||
          (nctype == NC_DOUBLE && mpitype == MPI_DOUBLE) );
#endif
}

static void
swapn(void *dst, const void *src, MPI_Offset nn, int xsize)
{
  int i;
  char *op = dst;
  const char *ip = src;
  while (nn-- != 0) {
    for (i=0; i<xsize; i++)
      op[i] = ip[xsize-1-i];
    op += xsize;
    ip += xsize;
  }
}


/* Endianness byte swap: done in-place */
#define SWAP(x,y) {tmp = (x); (x) = (y); (y) = tmp;}

/*----< in_swap() >----------------------------------------------------------*/
static void
in_swapn(void       *buf,
         MPI_Offset  nelems,  /* number of elements in buf[] */
         int         esize)   /* byte size of each element */
{
    int  i;
    char tmp, *op = buf;

    if (esize <= 1 || nelems <= 0) return;  /* no need */

    if (esize == 4) { /* this is the most common case */
        uint32_t *dest = (uint32_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = htonl(dest[i]);
    }
    else if (esize == 2) {
        uint16_t *dest = (uint16_t*) buf;
        for (i=0; i<nelems; i++)
            dest[i] = htons(dest[i]);
    }
    else {
        /* for esize is not 1, 2, or 4 */
        while (nelems-- > 0) {
            for (i=0; i<esize/2; i++)
                SWAP(op[i], op[esize-1-i])
            op += esize;
        }
    }
}


static int
x_putn_schar(void *xbuf, const void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int status = ENOERR;

  xp = (void *) xbuf;
  data = (void *) buf;


  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_SHORT)
      status = ncmpix_putn_schar_short(&xp, nelems, (const short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_putn_schar_int(&xp, nelems, (const int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_putn_schar_long(&xp, nelems, (const long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_putn_schar_float(&xp, nelems, (const float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_putn_schar_double(&xp, nelems, (const double *)data);

  return status;
}

static int
x_putn_short(void *xbuf, const void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int datainc;
  int status = ENOERR;
 
  xp = (void *) xbuf;
  data = (void *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_putn_short_uchar(&xp, nelems, (const unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_putn_short_schar(&xp, nelems, (const signed char *)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_putn_short_short(&xp, nelems, (const short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_putn_short_int(&xp, nelems, (const int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_putn_short_long(&xp, nelems, (const long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_putn_short_float(&xp, nelems, (const float *)data);
  else if (datatype == MPI_DOUBLE)
     status = ncmpix_putn_short_double(&xp, nelems, (const double *)data);

  return status;
} 

static int
x_putn_int(void *xbuf, const void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int datainc;
  int status = ENOERR;
 
  xp = (void *) xbuf;
  data = (void *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_putn_int_uchar(&xp, nelems, (const unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_putn_int_schar(&xp, nelems, (const signed char *)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_putn_int_short(&xp, nelems, (const short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_putn_int_int(&xp, nelems, (const int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_putn_int_long(&xp, nelems, (const long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_putn_int_float(&xp, nelems, (const float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_putn_int_double(&xp, nelems, (const double *)data);

  return status;
} 

static int
x_putn_float(void *xbuf, const void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int datainc;
  int status = ENOERR;
 
  xp = (void *) xbuf;
  data = (void *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_putn_float_uchar(&xp, nelems, (const unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_putn_float_schar(&xp, nelems, (const signed char *)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_putn_float_short(&xp, nelems, (const short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_putn_float_int(&xp, nelems, (const int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_putn_float_long(&xp, nelems, (const long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_putn_float_float(&xp, nelems, (const float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_putn_float_double(&xp, nelems, (const double *)data);

  return status;
} 

int
x_putn_double(void *xbuf, const void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int datainc;
  int status = ENOERR;
 
  xp = (void *) xbuf; 
  data = (void *) buf;
  datainc = length_of_mpitype(datatype);  

  if (datatype ==MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_putn_double_uchar(&xp, nelems, (const unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_putn_double_schar(&xp, nelems, (const signed char *)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_putn_double_short(&xp, nelems, (const short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_putn_double_int(&xp, nelems, (const int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_putn_double_long(&xp, nelems, (const long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_putn_double_float(&xp, nelems, (const float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_putn_double_double(&xp, nelems, (const double *)data);

  return status;
} 

static int
x_getn_schar(const void *xbuf, void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  void *xp, *data;
  int status = ENOERR;

  xp = (void *) xbuf;
  data = (void *) buf;

  if (datatype == MPI_CHAR)
        status = NC_ECHAR;
  else if (datatype == MPI_SHORT)
      status = ncmpix_getn_schar_short((const void **)&xp, nelems, (short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_getn_schar_int((const void **)&xp, nelems, (int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_getn_schar_long((const void **)&xp, nelems, (long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_getn_schar_float((const void **)&xp, nelems, (float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_getn_schar_double((const void **)&xp, nelems, (double *)data);

  return status;
}

static int
x_getn_short(const void *xbuf, void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
  int status = ENOERR;

  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_getn_short_uchar((const void **)&xp, nelems, (unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_getn_short_schar((const void **)&xp, nelems, (signed char*)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_getn_short_short((const void **)&xp, nelems, (short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_getn_short_int((const void **)&xp, nelems, (int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_getn_short_long((const void **)&xp, nelems, (long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_getn_short_float((const void **)&xp, nelems, (float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_getn_short_double((const void **)&xp, nelems, (double *)data);

  return status;
} 

static int 
x_getn_int(const void *xbuf, void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
  int status = ENOERR;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_getn_int_uchar((const void **)&xp, nelems, (unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_getn_int_schar((const void **)&xp, nelems, (signed char*)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_getn_int_short((const void **)&xp, nelems, (short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_getn_int_int((const void **)&xp, nelems, (int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_getn_int_long((const void **)&xp, nelems, (long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_getn_int_float((const void **)&xp, nelems, (float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_getn_int_double((const void **)&xp, nelems, (double *)data);

  return status;
} 

static int
x_getn_float(const void *xbuf, void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
  int  status = ENOERR;

  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_getn_float_uchar((const void **)&xp, nelems, (unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_getn_float_schar((const void **)&xp, nelems, (signed char*)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_getn_float_short((const void **)&xp, nelems, (short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_getn_float_int((const void **)&xp, nelems, (int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_getn_float_long((const void **)&xp, nelems, (long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_getn_float_float((const void **)&xp, nelems, (float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_getn_float_double((const void **)&xp, nelems, (double *)data);

  return status;
}

static int
x_getn_double(const void *xbuf, void *buf, MPI_Offset nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
  int  status = ENOERR;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
                                                                                    
  if (datatype == MPI_CHAR)
      status = NC_ECHAR;
  else if (datatype == MPI_UNSIGNED_CHAR)
      status = ncmpix_getn_double_uchar((const void **)&xp, nelems, (unsigned char *)data);
  else if (datatype == MPI_BYTE)
      status = ncmpix_getn_double_schar((const void **)&xp, nelems, (signed char*)data);
  else if (datatype == MPI_SHORT)
      status = ncmpix_getn_double_short((const void **)&xp, nelems, (short *)data);
  else if (datatype == MPI_INT)
      status = ncmpix_getn_double_int((const void **)&xp, nelems, (int *)data);
  else if (datatype == MPI_LONG)
      status = ncmpix_getn_double_long((const void **)&xp, nelems, (long *)data);
  else if (datatype == MPI_FLOAT)
      status = ncmpix_getn_double_float((const void **)&xp, nelems, (float *)data);
  else if (datatype == MPI_DOUBLE)
      status = ncmpix_getn_double_double((const void **)&xp, nelems, (double *)data);

  return status;
}

static int
check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm, int collective) {

  if (collective && NC_indep(ncp))
    return NC_EINDEP;

  if (!collective && !NC_indep(ncp))
    return NC_ENOTINDEP;

  if ( (collective && !NC_collectiveFhOpened(ncp->nciop)) 
    || (!collective && !NC_independentFhOpened(ncp->nciop)) ) {
  
    int mpireturn;
    mpireturn = MPI_File_open(comm, (char *)ncp->nciop->path, ncp->nciop->mpiomode,
                              ncp->nciop->mpiinfo, mpifh);
    if (mpireturn != MPI_SUCCESS) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      ncmpii_handle_error(rank, mpireturn, "check_mpifh(): MPI_File_open");
      return NC_ENFILE;
      /* To be determined the return error code ???????????? */
    }

    if (collective)
      set_NC_collectiveFh(ncp->nciop);
    else
      set_NC_independentFh(ncp->nciop);

  }

  return NC_NOERR;
}

static int check_recsize_too_big(NC *ncp) 
{
    int ret = NC_NOERR;
    /* assertion: because recsize will be used to set up the file
     * view, we must ensure there is no overflow when specifying
     * how big a stride there is between items (think interleaved
     * records).  
     *
     * note: 'recsize' is the sum of the record size of all record
     * variables in this dataset */
    if (ncp->recsize != (MPI_Aint)ncp->recsize) {
	    fprintf(stderr, "Type overflow: unable to read/write multiple records in this dataset\non this platform. Please either access records of this record variable\none-at-a-time or run on a 64 bit platform\n");
	    ret = NC_ESMALL;
    }
    /* the assert here might harsh, but without it, users will get corrupt
     * data.  */
    assert (ncp->recsize == (MPI_Aint)ncp->recsize);
    return ret;
}


/*----< NCcoordck() >--------------------------------------------------------*/
/*
 * Check whether 'coord' values (indices) are valid for the variable.
 */
static int
NCcoordck(NC               *ncp,
          const NC_var     *varp,
          const MPI_Offset *coord)
{
    const MPI_Offset *ip;
    MPI_Offset *up;
 
    if (varp->ndims == 0)
        return NC_NOERR;        /* 'scalar' variable */
 
    if (IS_RECVAR(varp)) {
/*      if (*coord > X_INT64_T_MAX)
            return NC_EINVALCOORDS; *//* sanity check */

        if (NC_readonly(ncp) && *coord >= ncp->numrecs) {
            if (!NC_doNsync(ncp))
                return NC_EINVALCOORDS;
            /* else */
            {
                /* Update from disk and check again */
                const int status = ncmpii_read_numrecs(ncp);
                if (status != NC_NOERR)
                    return status;
                if (*coord >= ncp->numrecs)
                    return NC_EINVALCOORDS;
            }
        }
        /* skip checking the record dimension */
        ip = coord + 1;
        up = varp->shape + 1;
    }
    else {
        ip = coord;
        up = varp->shape;
    }
 
    for (; ip < coord + varp->ndims; ip++, up++) {
        if ( (*ip <0) || (*ip >= *up) )
            return NC_EINVALCOORDS;
    }
    return NC_NOERR;
}

/*
 * Check whether 'edges' are valid for the variable and 'start'
 */
/*ARGSUSED*/
static int
NCedgeck(const NC *ncp, const NC_var *varp,
         const MPI_Offset *start, const MPI_Offset *edges)
{
  const MPI_Offset *const end = start + varp->ndims;
  const MPI_Offset *shp = varp->shape;

  if(varp->ndims == 0)
    return NC_NOERR;  /* 'scalar' variable */

  if(IS_RECVAR(varp))
  {
    start++;
    edges++;
    shp++;
  }

  for(; start < end; start++, edges++, shp++)
  {
    if( (*shp < 0) || (*edges > *shp) || (*start + *edges > *shp))
    {
      return(NC_EEDGE);
    }
  }

  return NC_NOERR;
}

#if 1 
/* enabled by Jianwei, 12/21/2004. why disable this? */
static int
NCstrideedgeck(const NC *ncp, const NC_var *varp,
         const MPI_Offset *start, const MPI_Offset *edges, const MPI_Offset *stride)
{
  const MPI_Offset *const end = start + varp->ndims;
  const MPI_Offset *shp = varp->shape; /* use MPI_Offset for now :( */

  if(varp->ndims == 0)
    return NC_NOERR;  /* 'scalar' variable */

  if(IS_RECVAR(varp))
  {
    if ( *stride == 0 ) /*|| *stride >= X_INT64_T_MAX)*/
      /* cast needed for braindead systems with signed MPI_Offset */
      return NC_ESTRIDE;

    start++;
    edges++;
    shp++;
    stride++;
  }

  for(; start < end; start++, edges++, shp++, stride++)
  {
    if( (*shp < 0) ||
        (*edges > *shp) || 
	(*edges > 0 && *start+1 + (*edges-1) * *stride > *shp) ||
	(*edges == 0 && *start > *shp) )
    {
      return(NC_EEDGE);
    }

    if ( *stride == 0)/* || *stride >= X_INT64_T_MAX)*/
      /* cast needed for braindead systems with signed MPI_Offset */
      return NC_ESTRIDE;
  }

  return NC_NOERR;
}
#endif

/*----< get_offset() >-------------------------------------------------------*/
static int
get_offset(NC               *ncp,
           NC_var           *varp,
           const MPI_Offset  starts[],   /* offsets relative to varp */
           MPI_Offset       *offset_ptr) /* return file offset */
{
    /* returns the starting file offset when this variable is get/put
       with starts[] */
    MPI_Offset offset;
    int status, dim, ndims;

    status = NCcoordck(ncp, varp, starts);  /* validate starts[] */
    if (status != NC_NOERR) {
        // printf("get_offset(): NCcoordck() fails\n");
        return status;
    }

    offset = varp->begin; /* beginning file offset of this variable */
    ndims  = varp->ndims; /* number of dimensions of this variable */

    if (ndims > 0) {
        if (IS_RECVAR(varp))
            /* no need to check recsize here: if MPI_Offset is only 32 bits we
               will have had problems long before here */
            offset += starts[0] * ncp->recsize;
        else
            offset += starts[ndims-1] * varp->xsz;

        if (ndims > 1) {
            if (IS_RECVAR(varp))
                offset += starts[ndims - 1] * varp->xsz;
            else
                offset += starts[0] * varp->dsizes[1] * varp->xsz;
    
            for (dim = 1; dim < ndims - 1; dim++)
                offset += starts[dim] * varp->dsizes[dim+1] * varp->xsz;
        }
    }
    *offset_ptr = offset;
    return NC_NOERR;
}

/*----< is_request_contiguous() >--------------------------------------------*/
static int
is_request_contiguous(NC_var           *varp,
                      const MPI_Offset  starts[],
                      const MPI_Offset  counts[])
{
    /* determine whether the get/put request to this variable using
       starts[] and counts[] is contiguous in file */
    int i, j, most_sig_dim, ndims = varp->ndims;

    /* this variable is a scalar */
    if (ndims == 0) return 1;

    most_sig_dim = 0; /* record dimension */
    if (IS_RECVAR(varp)) {
        /* if there are more than one record variabl, then the record
           dimensions, counts[0] must == 1. For now, we assume there
           are more than one record variable. 
           TODO: we need the API ncmpi_inq_rec() as in netcdf 3.6.3
                 to know how many record variables are defined */
        if (counts[0] != 1) return 0;
        most_sig_dim = 1;
    }

    for (i=ndims-1; i>most_sig_dim; i--) {
        /* find the first counts[i] that is not the entire dimension */
        if (counts[i] < varp->shape[i]) {
            /* check dim from i-1, i-2, ..., most_sig_dim and
               their counts[] should all be 1 */
            for (j=i-1; j>=most_sig_dim; j--) {
                if (counts[j] != 1) {
                    return 0;
                }
            }
            break;
        }
        else { /* counts[i] == varp->shape[i] */
            /* when accessing the entire dimension, starts[i] must be 0 */
            if (starts[i] != 0) return 0;
        }
    }
    return 1;
}

/*----< set_var_fileview() >-------------------------------------------------*/
static int
set_var_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp) {
    /* setting the fileview for entire variable. For non-record
       variable, the fileview is contiguous in file and hence
       no need the to create a filetype, but just use the starting
       offset */
    MPI_Offset offset;
    int mpireturn;

    offset = varp->begin;

    if (!IS_RECVAR(varp)) { 
        /* Contiguous file view */
        mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, MPI_BYTE,
                                      "native", MPI_INFO_NULL);
        if (mpireturn != MPI_SUCCESS) {
            int rank;
            MPI_Comm_rank(ncp->nciop->comm, &rank);
            ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
            return NC_EFILE;
        }
    }
    else {
        /* This is a record variable. If there are more than one record
           varaibles and this variable has more than one record, then the
           file view is strided (non-contiguous) in file */
        int  ndims;
        MPI_Datatype filetype;  
        MPI_Aint stride;
        int blocklen;

        if (ncp->numrecs == 0) /* no record been added yet */
            return(NC_NOERR);
        check_recsize_too_big(ncp);

        if (ncp->numrecs == 1) {
            /* there is only one record per record variable, the fileview
               is contiguous, then no need to create a filetype */
            mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE,
                                          MPI_BYTE, "native", MPI_INFO_NULL);
            if (mpireturn != MPI_SUCCESS) {
                int rank;
                MPI_Comm_rank(ncp->nciop->comm, &rank);
                ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
                return NC_EFILE;
            }
            return NC_NOERR;
        }

        /* now, this record variable has ncp->numrecs > 1 */
        ndims = varp->ndims;
        if (ndims > 1)
            blocklen = varp->dsizes[1] * varp->xsz;
        else
            blocklen = varp->xsz;

        /* ncp->recsize is the sum of 1st record of all record variables */
        stride = ncp->recsize;
    
#if (MPI_VERSION < 2)
        MPI_Type_hvector(ncp->numrecs, blocklen, stride, MPI_BYTE, &filetype);
#else
        MPI_Type_create_hvector(ncp->numrecs, blocklen, stride, MPI_BYTE, &filetype);
#endif
        MPI_Type_commit(&filetype);

        mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, filetype, "native", MPI_INFO_NULL);
        if (mpireturn != MPI_SUCCESS) {
            int rank;
            MPI_Comm_rank(ncp->nciop->comm, &rank);
            ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
            return NC_EFILE;
        }
        MPI_Type_free(&filetype); 
    }
    return NC_NOERR;
}

/*----< create_subarray_c_order_byte() >-------------------------------------*/
static int
create_subarray_c_order_byte(int ndims,
                             int sizes[],
                             int subsizes[],
                             int starts[],
                             MPI_Datatype *newtype)
{
    /* if a request is contiguous in file, this function calls
       MPI_Type_contiguous() and MPI_Type_create_resized(), we
       MPI_Type_create_subarray() can be avoided */
    int i, iscontig = 1;

    for (i=ndims-1; i>0; i--) {
        if (subsizes[i] < sizes[i]) {
            iscontig = 0;
            break;
        }
    }

    if (iscontig) {
        /* subsizes[1...ndims-1] == sizes[1...ndims-1] and
             starts[1...ndims-1] == 0 */
        MPI_Aint extent=sizes[0];

        *newtype = MPI_BYTE;
        for (i=ndims-1; i>0; i--) {
            MPI_Type_contiguous(sizes[i], *newtype, newtype);
            extent *= sizes[i];
        }
        MPI_Type_indexed(1, subsizes, starts, *newtype, newtype);
        /* augment the upper bound to the entire array size */
        return MPI_Type_create_resized(*newtype, 0, extent, newtype);
    }

    return MPI_Type_create_subarray(ndims, sizes, subsizes, starts,
                                    MPI_ORDER_C, MPI_BYTE, newtype);
}

/*----< set_vara_fileview() >------------------------------------------------*/
static int
set_vara_fileview(NC               *ncp,
                  MPI_File         *mpifh,
                  NC_var           *varp,
                  const MPI_Offset *start,
                  const MPI_Offset *count,
                  int               getnotput)
{
    int i, status, dim, ndims, mpireturn;
    int *shape = NULL, *subcount = NULL, *substart = NULL; /* all in bytes */
    MPI_Offset *shape64 = NULL, *subcount64 = NULL, *substart64 = NULL;
    MPI_Datatype rectype, filetype, types[3], type1;
    MPI_Offset offset, size, disps[3];
    int blklens[3], tag = 0;

    offset = varp->begin;
  
    ndims = varp->ndims;

    /* New coordinate/edge check to fix NC_EINVALCOORDS bug */
    status = NCedgeck(ncp, varp, start, count);
    if (status != NC_NOERR ||
        (getnotput && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp)) ) 
    {
        status = NCcoordck(ncp, varp, start);
        if (status != NC_NOERR)
            return status;
        else
            return NC_EEDGE;
    }

/* Removed to fix NC_EINVALCOORDS bug

    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR)
        return status;

    status = NCedgeck(ncp, varp, start, count);
    if (status != NC_NOERR)
        return status;

    if (getnotput && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp))
        return NC_EEDGE;
*/

    /* check if the request is contiguous in file
       if yes, there is no need to create a filetype */
    if (is_request_contiguous(varp, start, count)) {
        status = get_offset(ncp, varp, start, &offset);
        if (status != NC_NOERR) return status;

        mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, 
                                      MPI_BYTE, "native", MPI_INFO_NULL);
        if (mpireturn != MPI_SUCCESS) {
            int rank;
            MPI_Comm_rank(ncp->nciop->comm, &rank);
            ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
            return NC_EFILE;
        }
        return NC_NOERR;
    }

    if (ndims == 0) {
        /* scalar variable */
        filetype = MPI_BYTE;
    }
    else {
        /* if ndims == 0, all below pointers would be null */
        shape    = (int *)malloc(sizeof(int) * ndims);
        subcount = (int *)malloc(sizeof(int) * ndims);
        substart = (int *)malloc(sizeof(int) * ndims);

        dim = 0;
        while (dim < ndims && count[dim] > 0) dim++;

        if (dim < ndims) {
            /* 0 size data */
            filetype = MPI_BYTE;
        }
        else {
           if (IS_RECVAR(varp)) {
               subcount[0] = count[0];
               substart[0] = 0;
               shape[0] = subcount[0];

               if (ncp->recsize <= varp->len) {
                   /* the only record variable */

                   if (varp->ndims == 1) {
                       shape[0] *= varp->xsz;
                       subcount[0] *= varp->xsz;
                   }
                   else {
                       for (dim = 1; dim < ndims-1; dim++) {
                           shape[dim]    = varp->shape[dim];
                           subcount[dim] = count[dim];
                           substart[dim] = start[dim];
                       }
                       shape[dim]    = varp->xsz * varp->shape[dim];
                       subcount[dim] = varp->xsz * count[dim];
                       substart[dim] = varp->xsz * start[dim];
                   }
                   offset += start[0] * ncp->recsize;

                   MPI_Type_create_subarray(ndims, shape, subcount, substart, 
                                            MPI_ORDER_C, MPI_BYTE, &filetype); 
                   MPI_Type_commit(&filetype);
               }
               else {
                   check_recsize_too_big(ncp);
                   /* more than one record variables */

                   offset += start[0] * ncp->recsize;
                   if (varp->ndims == 1) {
#if (MPI_VERSION < 2)
                       MPI_Type_hvector(subcount[0], varp->xsz, ncp->recsize,
                                        MPI_BYTE, &filetype);
#else
                       MPI_Type_create_hvector(subcount[0], varp->xsz, ncp->recsize,
                                               MPI_BYTE, &filetype);
#endif
                       MPI_Type_commit(&filetype);
                   }
                   else {
                       for (dim = 1; dim < ndims-1; dim++) {
                           shape[dim]    = varp->shape[dim];
                           subcount[dim] = count[dim];
                           substart[dim] = start[dim];
                       }
                       shape[dim]    = varp->xsz * varp->shape[dim];
                       subcount[dim] = varp->xsz * count[dim];
                       substart[dim] = varp->xsz * start[dim];

                       MPI_Type_create_subarray(ndims-1, shape+1, subcount+1, substart+1,
                                                MPI_ORDER_C, MPI_BYTE, &rectype);
                       MPI_Type_commit(&rectype);
#if (MPI_VERSION < 2)
                       MPI_Type_hvector(subcount[0], 1, ncp->recsize, rectype, &filetype);
#else
                       MPI_Type_create_hvector(subcount[0], 1, ncp->recsize, rectype, &filetype);
#endif
                       MPI_Type_commit(&filetype);
                       MPI_Type_free(&rectype);
                   }
               }
           }
           else { /* non record variable */
               tag = 0;
               for (dim=0; dim< ndims-1; dim++) {
                   if (varp->shape[dim] > 2147483647) { /* if shape > 2^31-1 */
                       tag = 1;
                        break;
                   }
               }
               if ((varp->shape[dim]*varp->xsz)  > 2147483647)
                   tag = 1;

               if (tag == 0) {
                   for (dim = 0; dim < ndims-1; dim++ ) {
                       shape[dim]    = varp->shape[dim];
                       subcount[dim] = count[dim];
                       substart[dim] = start[dim];
                   } 
                   shape[dim]    = varp->xsz * varp->shape[dim];
                   subcount[dim] = varp->xsz * count[dim];
                   substart[dim] = varp->xsz * start[dim];
  
                   MPI_Type_create_subarray(ndims, shape, subcount, substart, 
                                            MPI_ORDER_C, MPI_BYTE, &filetype); 
                   MPI_Type_commit(&filetype);
               }
               else {
                   shape64    = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims);
                   subcount64 = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims);
                   substart64 = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims);

                   if (ndims == 1) {  // for 64-bit support,  added July 23, 2008
                       shape64[0]    = varp->shape[0];
                       subcount64[0] = count[0];
                       substart64[0] = start[0];

                       offset += start[0]*varp->xsz;

                       MPI_Type_contiguous(subcount64[0]*varp->xsz, MPI_BYTE, &type1);
                       MPI_Type_commit(&type1);
#if (MPI_VERSION < 2)
                       MPI_Type_hvector(subcount64[0], varp->xsz, shape64[0]*varp->xsz,
                                        MPI_BYTE, &filetype);
#else
                       MPI_Type_create_hvector(1, 1, shape64[0]*varp->xsz,
                                               type1, &filetype);
#endif
                       MPI_Type_commit(&filetype);
                       MPI_Type_free(&type1);
                   }
                   else {
                       for (dim = 0; dim < ndims-1; dim++ ) {
                           shape64[dim]    = varp->shape[dim];
                           subcount64[dim] = count[dim];
                           substart64[dim] = start[dim];
                       }
                       shape64[dim]    = varp->xsz * varp->shape[dim];
                       subcount64[dim] = varp->xsz * count[dim];
                       substart64[dim] = varp->xsz * start[dim];

                       MPI_Type_hvector(subcount64[dim-1],
                                        subcount64[dim],
                                        varp->xsz * varp->shape[dim],
                                        MPI_BYTE,
                                        &type1);
                       MPI_Type_commit(&type1);

                       size = shape[dim];
                       for (i=dim-2; i>=0; i--) {
                           size *= shape[i+1];
                           MPI_Type_hvector(subcount64[i],
                                            1,
                                            size,
                                            type1,
                                            &filetype);
                           MPI_Type_commit(&filetype);

                           MPI_Type_free(&type1);
                           type1 = filetype;
                       }
                       disps[1] = substart64[dim];
                       size = 1;
                       for (i=dim-1; i>=0; i--) {
                           size *= shape64[i+1];
                           disps[1] += size*substart64[i];
                       }
                       disps[2] = 1;
                       for (i=0; i<ndims; i++) disps[2] *= shape64[i];

                       disps[0] = 0;
                       blklens[0] = blklens[1] = blklens[2] = 1;
                       types[0] = MPI_LB;
                       types[1] = type1;
                       types[2] = MPI_UB;

                       MPI_Type_struct(3,
                                       blklens,
                                       (MPI_Aint*) disps,
                                       types,
                                       &filetype);

                       MPI_Type_free(&type1);
                   }
                   free(shape64);
                   free(subcount64);
                   free(substart64);
               }
           }
       }
   }

   mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, 
                                 filetype, "native", MPI_INFO_NULL);
   if (mpireturn != MPI_SUCCESS) {
       int rank;
       MPI_Comm_rank(ncp->nciop->comm, &rank);
       ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
       return NC_EFILE;
   }

   if (ndims > 0) {
       if (filetype != MPI_BYTE)
           MPI_Type_free(&filetype);
           
       free(shape);
       free(subcount);
       free(substart);
   }
   return NC_NOERR;
}

static int
set_vars_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, 
		     const MPI_Offset start[], const MPI_Offset count[], 
                     const MPI_Offset stride[], int getnotput) {
  MPI_Offset offset;
  int status;
  int mpireturn;
  int dim, ndims;
  MPI_Datatype *subtypes, *filetype;
  MPI_Offset *blocklens = NULL, *blockstride = NULL, *blockcount = NULL;

  ndims = varp->ndims;

  for (dim=0; dim<ndims && stride[dim]==1; dim++) ;
  if (dim == ndims)
    return set_vara_fileview(ncp, mpifh, varp, start, count, getnotput);

  /* New coordinate/edge check to fix NC_EINVALCOORDS bug */
  status = NCedgeck(ncp, varp, start, count);
  if( (status != NC_NOERR) ||
      (getnotput && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp)) )
  {
    status = NCcoordck(ncp, varp, start);
    if (status != NC_NOERR)
      return status;
    else
      return NC_EEDGE;
  }

  status = NCstrideedgeck(ncp, varp, start, count, stride);
  if(status != NC_NOERR)
    return status;

  if( getnotput && IS_RECVAR(varp) && 
     ( (*count > 0 && *start+1 + (*count-1) * *stride > NC_get_numrecs(ncp)) ||
       (*count == 0 && *start > NC_get_numrecs(ncp)) ) )
    return NC_EEDGE;

/* Removed to fix NC_EINVALCOORDS bug 

  status = NCcoordck(ncp, varp, start);
  if (status != NC_NOERR)
    return status; 
*/

/* Moved into NCstrideedgeck

  for (dim = 0; dim < ndims; dim++)
  {
    if ( (stride != NULL && stride[dim] == 0) ||
        stride[dim] >= X_LONG_MAX)
    {
      return NC_ESTRIDE;
    }
  }
*/

/* Removed to fix NC_EINVALCOORDS bug 

  status = NCedgeck(ncp, varp, start, count);
  if(status != NC_NOERR)
    return status;

 if(getnotput && IS_RECVAR(varp) &&
     (MPI_Offset)*start + (MPI_Offset)*count > NC_get_numrecs(ncp))
      return NC_EEDGE;
*/

  offset = varp->begin;
  
  if (ndims == 0) {

    /* scalar variable */

    filetype = subtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype));
    *filetype = MPI_BYTE;
  
  } else {  

    blocklens = (MPI_Offset *) malloc(ndims * sizeof(MPI_Offset));
    blockstride = (MPI_Offset *) malloc(ndims * sizeof(MPI_Offset));
    blockcount = (MPI_Offset *) malloc(ndims * sizeof(MPI_Offset));

    dim = 0;
    while (dim < ndims && count[dim] > 0) dim++;

    if (dim < ndims) {

      /* 0 size data */
      filetype = subtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype));
      *filetype = MPI_BYTE;

    } else {
  
      subtypes = (MPI_Datatype *)malloc((ndims+1) * sizeof(MPI_Datatype));
      filetype = subtypes;
      subtypes[ndims] = MPI_BYTE;
  
      blocklens[ndims - 1] = varp->xsz;
      blockcount[ndims - 1] = count[ndims - 1];
      if (ndims == 1 && IS_RECVAR(varp)) {
	check_recsize_too_big(ncp);
        blockstride[ndims - 1] = stride[ndims - 1] * ncp->recsize;
        offset += start[ndims - 1] * ncp->recsize;
      } else {
        blockstride[ndims - 1] = stride[ndims - 1] * varp->xsz;
        offset += start[ndims - 1] * varp->xsz;
      }
  
      for (dim = ndims - 1; dim >= 0; dim--) {
#if (MPI_VERSION < 2)
        MPI_Type_hvector(blockcount[dim], blocklens[dim], blockstride[dim],
                         subtypes[dim + 1], subtypes + dim);
#else
        MPI_Type_create_hvector(blockcount[dim], blocklens[dim], blockstride[dim],
                              subtypes[dim + 1], subtypes + dim);
#endif
        MPI_Type_commit(subtypes + dim);
      
        if (dim - 1 >= 0) {
          blocklens[dim - 1] = 1;
          blockcount[dim - 1] = count[dim - 1];
          if (dim - 1 == 0 && IS_RECVAR(varp)) {
            blockstride[dim - 1] = stride[dim - 1] * ncp->recsize;
	    offset += start[dim-1] * ncp->recsize;
          } else {
            blockstride[dim - 1] = stride[dim - 1] * varp->dsizes[dim] * varp->xsz;
            offset += start[dim-1] * varp->dsizes[dim] * varp->xsz;
          }
        }
      } 

    }

  }

  mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, 
                    *filetype, "native", MPI_INFO_NULL);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
        return NC_EFILE;
  }

  if (ndims > 0) {
    if (*filetype != MPI_BYTE)
      for (dim=0; dim < ndims; dim++) 
        MPI_Type_free(subtypes + dim);

    free(blocklens);
    free(blockstride);
    free(blockcount);
  }

  free(subtypes);

  return NC_NOERR;
}

/* BEGIN put/get functions */
/* buffer layers:	
	
	User Level		buf	(user defined buffer of MPI_Datatype)
	MPI Datatype Level	cbuf	(contiguous buffer of ptype)
	NetCDF XDR Level	xbuf	(XDR I/O buffer)
*/

int
ncmpi_put_var1(int ncid, int varid,
               const MPI_Offset index[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes, offset;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;

  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status; 

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  if (nelems != cnelems) {
    if (warning == NC_NOERR) 
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems*varp->xsz; /* account for file bytes */
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* accessing one variable element need not set the file view.
     just find the file offset and use MPI-IO call with explicit offset */

  status = get_offset(ncp, varp, index, &offset);
  if (status != NC_NOERR) return status;

  if (!iscontig_of_ptypes) {
  
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype, 
				cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      goto fn_exit;
 
  } else {
   
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */
 
  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */
  
    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, 1, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, 1, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, 1, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, 1, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, 1, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, 1, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write_at(ncp->nciop->independent_fh, offset, xbuf, nbytes,
			     MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write_at");
        status = NC_EWRITE;
  }
 
  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
      in_swapn(cbuf, 1, ncmpix_len_nctype(varp->type));

fn_exit:
  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if ( status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = index[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_var1(int ncid, int varid,
               const MPI_Offset index[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes, offset; 
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status; 

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems*varp->xsz; /* account for file bytes */
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* it is not necessary to define a filetype because requesting one
     array element is accessing a contiguous region in file */
  status = get_offset(ncp, varp, index, &offset);
  if (status != NC_NOERR) return status;

  if (!iscontig_of_ptypes) {
  
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read_at(ncp->nciop->independent_fh, offset, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read_at");
        status = NC_EREAD;
  }
 
  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, 1, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, 1, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, 1, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, 1, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, 1, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, 1, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype, 
				(void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_var_all(int ncid, int varid, void *buf, MPI_Offset bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  if (varp->ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (varp->ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems * varp->xsz;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* since reading an entire non-record variable is a contiguous file
     access, no need to set fileview. Setting fileview is only necessary
     for record variable and only when there are more than one record */
  if (IS_RECVAR(varp)) {
      /* Record variable has a strided file view */
      status = set_var_fileview(ncp, &(ncp->nciop->collective_fh), varp);
      if (status != NC_NOERR)
          return status;
  }
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }


  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  if (IS_RECVAR(varp))
      mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes,
                                    MPI_BYTE, &mpistatus);
  else
      /* Contiguous file view - read the entire variable from an offset */
      mpireturn = MPI_File_read_at_all(ncp->nciop->collective_fh, varp->begin,
                                       xbuf, nbytes, MPI_BYTE, &mpistatus);

  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read(_at)_all");
        status = NC_EREAD;
  }

  if (IS_RECVAR(varp)) /* reset fileview so the entire file is visible again */
      MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                        "native", MPI_INFO_NULL);
 
  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_var(int ncid, int varid, const void *buf, MPI_Offset bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  if (varp->ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (varp->ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems * varp->xsz;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* since writing an entire non-record variable is a contiguous file
     access, no need to set fileview. Setting fileview is only necessary
     for record variable and only when there are more than one record */
  if (IS_RECVAR(varp)) {
      /* Record variable has a strided file view */
      status = set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
      if (status != NC_NOERR)
          return status;
  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  if (IS_RECVAR(varp))
      mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes,
                                 MPI_BYTE, &mpistatus);
  else
      /* Contiguous file view - write the entire variable from an offset */
      mpireturn = MPI_File_write_at(ncp->nciop->independent_fh, varp->begin,
                                    xbuf, nbytes, MPI_BYTE, &mpistatus);

  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write(_at)");
        status = NC_EWRITE;
  }

  if (IS_RECVAR(varp)) /* reset fileview so the entire file is visible again */
      MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                        "native", MPI_INFO_NULL);
 
  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    if (varp->ndims > 1)
      newnumrecs = nelems / varp->dsizes[1];
    else
      newnumrecs = nelems;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_var(int ncid, int varid, void *buf, MPI_Offset bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  if (varp->ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (varp->ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems * varp->xsz;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  if (IS_RECVAR(varp)) {
      /* Accessing the entire variable, only record variables have a
         strided file view */
      status = set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
      if (status != NC_NOERR)
          return status;
  }

  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  if (IS_RECVAR(varp))
      mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes,
                                MPI_BYTE, &mpistatus);
  else
      /* Contiguous file view - read the entire variable from an offset */
      mpireturn = MPI_File_read_at(ncp->nciop->independent_fh, varp->begin,
                                   xbuf, nbytes, MPI_BYTE, &mpistatus);

  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read(_at)");
        status = NC_EREAD;
  }

  if (IS_RECVAR(varp)) /* reset fileview so the entire file is visible again */
      MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                        "native", MPI_INFO_NULL);
 
  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_put_vara_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount, 
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* set the mpi file view */
 
  status = set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp, start, count, 0);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {

    /* handling for derived datatype: pack into a contiguous buffer */

    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype, 
				cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
  
  } else {

    cbuf = (void *)buf;

  }


  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  /* if record variables are too big (so big that we cannot store the stride
   * between records in an MPI_Aint, for example) then we will have to process
   * this one record at a time.  
   *
   * It stinks that we have to make this change in multiple places by the way
   */ 

  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes,
                                 MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write_all");
        status = NC_EWRITE;
  }

  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);

  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
      /* update the number of records in header NC */
      MPI_Offset newnumrecs = start[0] + count[0];
      update_numrecs(ncp, newnumrecs);
  }
  return ((warning != NC_NOERR) ? warning : status);
}


int
ncmpi_get_vara_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
		   void *buf, MPI_Offset bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems; 
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* set the mpi file view */
 
  status = set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp, start, count, 1);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes,
                                MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read_all");
        status = NC_EREAD;
  }

  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);

  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion and byte swap if needed */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(xbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype, 
				(void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

/* encode 'buf' into netcdf format, both in terms of type (the 'type') and
 * endianess */
int ncmpiii_data_encode(MPI_Offset cnelems, int el_size, 
        int iscontig_of_ptypes, 
        const void *buf, void **cbufp, void **xbufp, 
        MPI_Offset bufcount, MPI_Datatype datatype,
        MPI_Datatype ptype,
        nc_type type, MPI_Offset nbytes, MPI_Offset nelems)
{ 
    int status;

    void *cbuf, *xbuf;
    cbuf = *cbufp;
    xbuf = *xbufp;

    if (!iscontig_of_ptypes) {

        /* handling for derived datatype: pack into a contiguous buffer */

        cbuf = (void *)malloc( cnelems * el_size );
        status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                cbuf, cnelems, ptype);
        if (status != NC_NOERR)
            return status;

    } else {

        cbuf = (void *)buf;

    }

    /* assign or allocate MPI buffer */

    if ( need_convert(type, ptype) ) {

        /* allocate new buffer */

        xbuf = (void *)malloc(nbytes);

        /* automatic numeric datatype conversion */

        switch( type ) {
            case NC_BYTE:
                status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_SHORT:
                status = x_putn_short(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_INT:
                status = x_putn_int(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_FLOAT:
                status = x_putn_float(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_DOUBLE:
                status = x_putn_double(xbuf, cbuf, cnelems, ptype);
                break;
            default:
                break;
        }

    } else if ( need_swap(type, ptype) ) {

        /* allocate new buffer */
        xbuf = (void *)malloc(nbytes);

        swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(type));

    } else {

        /* else, just assign contiguous buffer */
        xbuf = (void *)cbuf;

    }

    *xbufp = xbuf;
    *cbufp = cbuf;
    return NC_NOERR;
}

/* slightly different from ncmpiii_data_encode in that it operates on an
 * already-allocated buffer. it would be nice to fix that so these two routines
 * behave the same...  */

/* decode data from netcdf format in xbuf into user's format in 'buf' */
int ncmpiii_data_decode(MPI_Offset cnelems, int el_size, 
        int iscontig_of_ptypes, 
        void *buf, void *cbuf, void *xbuf, 
        MPI_Offset bufcount, MPI_Datatype datatype,
        MPI_Datatype ptype,
        nc_type type, MPI_Offset nbytes, MPI_Offset nelems)
{ 
    int status;


    if ( need_convert(type, ptype) ) {

        /* automatic numeric datatype conversion */

        switch( type ) {
            case NC_BYTE:
                status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_SHORT:
                status = x_getn_short(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_INT:
                status = x_getn_int(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_FLOAT:
                status = x_getn_float(xbuf, cbuf, cnelems, ptype);
                break;
            case NC_DOUBLE:
                status = x_getn_double(xbuf, cbuf, cnelems, ptype);
                break;
            default:
                break;
        }

    } else if ( need_swap(type, ptype) ) {

        if (xbuf == cbuf)
            in_swapn(cbuf, nelems, ncmpix_len_nctype(type));
        else
            swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(type));
    }
    if (!iscontig_of_ptypes) {
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype, 
				(void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
    }
 
    return NC_NOERR;
}

int
ncmpi_put_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL,  *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes, offset;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* check if the request is contiguous in file
     if yes, we will later use offset in MPI_File_write_at() */
  offset = -1;
  if (is_request_contiguous(varp, start, count)) {
      status = get_offset(ncp, varp, start, &offset);
      if (status != NC_NOERR) return status;
  }
  else {
      /* this request is non-contiguous in file, set the mpi file view */
      status = set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp,
                                 start, count, 0);
      if(status != NC_NOERR)
          return status;
  }
 
  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  if (offset >= 0) /* this is a contiguous file access */
      mpireturn = MPI_File_write_at(ncp->nciop->independent_fh, offset, xbuf,
                                    nbytes, MPI_BYTE, &mpistatus);
  else
      mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes,
                                 MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write(_at)");
        return NC_EWRITE;
  }

  if (offset < 0)
      /* reset the file view so the entire file is visible again */
      MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                        "native", MPI_INFO_NULL);

  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = start[0] + count[0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_get_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes, offset;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* check if the request is contiguous in file
     if yes, we will later use offset in MPI_File_read_at() */
  offset = -1;
  if (is_request_contiguous(varp, start, count)) {
      status = get_offset(ncp, varp, start, &offset);
      if (status != NC_NOERR) return status;
  }
  else {
      /* this request is non-contiguous in file, set the mpi file view */
      status = set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp,
                                 start, count, 1);
      if (status != NC_NOERR)
          return status;
  }
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }
 
  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  if (offset >= 0) /* this is a contiguous file access */
      mpireturn = MPI_File_read_at(ncp->nciop->independent_fh, offset, xbuf,
                                   nbytes, MPI_BYTE, &mpistatus);
  else
      mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes,
                                MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read(_at)");
        status = NC_EREAD;
  }

  if (offset < 0)
      /* reset the file view so the entire file is visible again */
      MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                        "native", MPI_INFO_NULL);

  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_put_vars_all(int ncid, int varid,
                   const MPI_Offset start[], 
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
                   const void *buf, MPI_Offset bufcount, 
                   MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
  MPI_Offset *stride_was_null=NULL, *stride_ptr;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* NULL stride is legal and means (1,1,1,..1) */
  stride_ptr = (MPI_Offset *)stride;
  if (stride == NULL) {
	  stride_was_null=(MPI_Offset*)malloc(varp->ndims*sizeof(MPI_Offset));
	  for (dim=0; dim < varp->ndims; dim++) {
		  stride_was_null[dim] = 1;
	  }
	  stride_ptr = stride_was_null;
  }

  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->collective_fh), 
                                varp, start, count, stride_ptr, 0);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes,
                                 MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write_all");
        status = NC_EWRITE;
  }

  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);

  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
      /* update the number of records in header NC */
      MPI_Offset newnumrecs = start[0] + (count[0] - 1) * stride_ptr[0] + 1;
      update_numrecs(ncp, newnumrecs);
  }
  if (stride_was_null != NULL) free(stride_was_null);
 
  return ((warning != NC_NOERR) ? warning : status);
}


int
ncmpi_get_vars_all(int ncid, int varid,
                   const MPI_Offset start[], 
		   const MPI_Offset count[],
                   const MPI_Offset stride[],
		   void *buf, MPI_Offset bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems; 
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->collective_fh), 
				varp, start, count, stride, 1);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes,
                                MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read_all");
        status = NC_EREAD;
  }

  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);

  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_vars(int ncid, int varid,
               const MPI_Offset start[], 
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
			        varp, start, count, stride, 0);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, datatype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes,
                             MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write");
        status = NC_EWRITE;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);

  if ( need_swap(varp->type, ptype) && cbuf == buf && cbuf == xbuf )
    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = start[0] + (count[0] - 1) * stride[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_get_vars(int ncid, int varid,
               const MPI_Offset start[], 
	       const MPI_Offset count[],
               const MPI_Offset stride[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
				varp, start, count, stride, 1); 
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes,
                            MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read");
        status = NC_EREAD;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE,
                    "native", MPI_INFO_NULL);
 
  if ( need_convert(varp->type, ptype) ) {

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));

  }

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: unpack from the contiguous buffer */
 
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;
 
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
} 

/* varm: there maybe two layer of memory layout (remapping):
	 one is specified by MPI derived datatype,
	 the other is specified by imap[],
   	 it's encouraged to use only one option of them,
	 though using both of them are supported.

   user buffer:				|--------------------------|

   mpi derived datatype view:		|------|  |------|  |------|
		
   logic (contig) memory datastream:	   |------|------|------|

   imap view:				   |--| |--|    |--| |--|

   contig I/O datastream (internal represent): |--|--|--|--|

   These two layers of memory layout will both be represented in MPI 
   derived datatype, and if double layers of memory layout is used, 
   we need to elimilate the upper one passed in MPI_Datatype parameter
   from the user, by repacking it to logic contig memory datastream view.
*/

int
ncmpi_put_varm_all(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   const void *buf, MPI_Offset bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems;
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                              buf, bufcount, datatype);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                              buf, bufcount, datatype);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                              buf, bufcount, datatype);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {

    /* handling for derived datatype: pack into a contiguous buffer */

    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      goto fn_exit;

  } else {

    lbuf = (void *)buf;

  }

  if (count[dim] < 0) {
    status = NC_ENEGATIVECNT;
    goto fn_exit;
  }
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen*count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0) {
      status = NC_ENEGATIVECNT;
      goto fn_exit;
    }

#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, imap[dim]*el_size,
                            imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  /* layout lbuf to cbuf based on imap */
  status = ncmpii_data_repack(lbuf, 1, imaptype,
			      cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    goto fn_exit;

  MPI_Type_free(&imaptype);

  status = ncmpi_put_vars_all(ncid, varid, start, count, stride,
                              cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    goto fn_exit;

fn_exit:
  if (!iscontig_of_ptypes && lbuf != NULL)
    free(lbuf);
  if (cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_varm_all(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   void *buf, MPI_Offset bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems; 
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                              buf, bufcount, datatype);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_get_vars_all(ncid, varid, start, count, stride,
			      buf, bufcount, datatype);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_get_vars_all(ncid, varid, start, count, stride,
			      buf, bufcount, datatype);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      goto fn_exit;
 
  } else {
 
    lbuf = (void *)buf;
 
  }

  if (count[dim] < 0) {
    status = NC_ENEGATIVECNT;
    goto fn_exit;
  }
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
		  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen * count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0) {
      status = NC_ENEGATIVECNT;
      goto fn_exit;
    }

#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, (MPI_Aint)imap[dim]*el_size, 
			    imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  status = ncmpi_get_vars_all(ncid, varid, start, count, stride, 
			      cbuf, cnelems, ptype);
  if (status != NC_NOERR) {
    if (status == NC_ERANGE && warning == NC_NOERR) 
      warning = status; /* to satisfy the nc_test logic */
    else
      goto fn_exit;

  }

  /* layout cbuf to lbuf based on imap */
  status = ncmpii_data_repack(cbuf, cnelems, ptype,
			      lbuf, 1, imaptype);
  if (status != NC_NOERR)
    goto fn_exit;

  MPI_Type_free(&imaptype);

  if (!iscontig_of_ptypes) {

    /* repack it back, like a read-modify-write operation */

    status = ncmpii_data_repack(lbuf, lnelems, ptype,
				(void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      goto fn_exit;

  }

fn_exit:
  if (!iscontig_of_ptypes && lbuf != NULL)
      free(lbuf);
  if (cbuf != NULL)
    free(cbuf);
  if (imaptype != MPI_DATATYPE_NULL)
	  MPI_Type_free(&imaptype);
    
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_varm(int ncid, int varid,
	       const MPI_Offset start[],
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
	       const MPI_Offset imap[],
	       const void *buf, MPI_Offset bufcount,
	       MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems;
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_put_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_put_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_put_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {

    /* handling for derived datatype: pack into a contiguous buffer */

    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      goto fn_exit;

  } else {

    lbuf = (void *)buf;

  }

  if (count[dim] < 0) {
    status = NC_ENEGATIVECNT;
    goto fn_exit;
  }
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen*count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0) {
      status = NC_ENEGATIVECNT;
      goto fn_exit;
    }
  
#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, imap[dim]*el_size,
                            imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  /* layout lbuf to cbuf based on imap */
  status = ncmpii_data_repack(lbuf, 1, imaptype,
			      cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    goto fn_exit;

  MPI_Type_free(&imaptype);

  status = ncmpi_put_vars(ncid, varid, start, count, stride,
                          cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    goto fn_exit;

fn_exit:
  if (!iscontig_of_ptypes && lbuf != NULL)
    free(lbuf);
  if (cbuf != NULL)
    free(cbuf);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_varm(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   void *buf, MPI_Offset bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems;
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_get_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_get_vars(ncid, varid, start, count, stride,
			  buf, bufcount, datatype);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_get_vars(ncid, varid, start, count, stride,
			  buf, bufcount, datatype);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    lbuf = (void *)buf;
 
  }

  if (count[dim] < 0)
    return NC_ENEGATIVECNT;
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
		  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen * count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
  
#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, (MPI_Aint)imap[dim]*el_size, 
			    imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  status = ncmpi_get_vars(ncid, varid, start, count, stride, 
			  cbuf, cnelems, ptype);
  if (status != NC_NOERR) {
    if (status == NC_ERANGE && warning == NC_NOERR) 
      warning = status; /* to satisfy the nc_test logic */
    else
      return status;
  }

  /* layout cbuf to lbuf based on imap */
  status = ncmpii_data_repack(cbuf, cnelems, ptype,
			      lbuf, 1, imaptype);
  if (status != NC_NOERR)
    return status;

  MPI_Type_free(&imaptype);

  if (!iscontig_of_ptypes) {

    /* repack it back, like a read-modify-write operation */

    status = ncmpii_data_repack(lbuf, lnelems, ptype,
				(void *)buf, bufcount, datatype);
    if (status != NC_NOERR)
      return status;

    if (lbuf != NULL)
      free(lbuf);

  }

  if (cbuf != NULL)
    free(cbuf);
    
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_var1_uchar(int ncid, int varid,
                     const MPI_Offset index[],
                     const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_var1_schar(int ncid, int varid,
                     const MPI_Offset index[],
                     const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_BYTE);
}

int
ncmpi_put_var1_text(int ncid, int varid,
                     const MPI_Offset index[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_CHAR);
}


int
ncmpi_put_var1_short(int ncid, int varid,
                     const MPI_Offset index[],
		     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index, 
                        (const void *)op, 1, MPI_SHORT); 
}

int
ncmpi_put_var1_int(int ncid, int varid,
                   const MPI_Offset index[],
                   const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_INT);
}

int
ncmpi_put_var1_long(int ncid, int varid,
                   const MPI_Offset index[],
                   const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_LONG);
}

int
ncmpi_put_var1_float(int ncid, int varid,
                     const MPI_Offset index[],
                     const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_FLOAT); 
}
 
int
ncmpi_put_var1_double(int ncid, int varid,
                      const MPI_Offset index[],
                      const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_DOUBLE);
}

int
ncmpi_get_var1_uchar(int ncid, int varid,
                     const MPI_Offset index[],
                     unsigned char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_var1_schar(int ncid, int varid,
                     const MPI_Offset index[],
                     signed char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_BYTE);
}

int
ncmpi_get_var1_text(int ncid, int varid,
                     const MPI_Offset index[],
                     char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_CHAR);
}

int
ncmpi_get_var1_short(int ncid, int varid,
                     const MPI_Offset index[],
                     short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_SHORT); 
}
 
int
ncmpi_get_var1_int(int ncid, int varid,
                   const MPI_Offset index[],
                   int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_INT);  
}
 
int
ncmpi_get_var1_long(int ncid, int varid,
                   const MPI_Offset index[],
                   long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_LONG); 
}

int
ncmpi_get_var1_float(int ncid, int varid,
                     const MPI_Offset index[],
                     float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_FLOAT);  
}
 
int
ncmpi_get_var1_double(int ncid, int varid,
                      const MPI_Offset index[],
                      double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_DOUBLE);  
} 

int
ncmpi_put_var_uchar(int ncid, int varid, const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_var_schar(int ncid, int varid, const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_BYTE);
}


int
ncmpi_put_var_text(int ncid, int varid, const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */
  
  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;
  
  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_var_short(int ncid, int varid, const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_var_int(int ncid, int varid, const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_INT);
} 

int
ncmpi_put_var_long(int ncid, int varid, const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_var_float(int ncid, int varid, const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_FLOAT);
} 

int
ncmpi_put_var_double(int ncid, int varid, const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nelems, MPI_DOUBLE);
} 

int
ncmpi_get_var_uchar(int ncid, int varid, unsigned char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_var_schar(int ncid, int varid, signed char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_var_text(int ncid, int varid, char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_var_short(int ncid, int varid, short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_var_int(int ncid, int varid, int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_INT);
} 

int
ncmpi_get_var_long(int ncid, int varid, long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_var_float(int ncid, int varid, float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_FLOAT);
} 

int
ncmpi_get_var_double(int ncid, int varid, double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nelems, MPI_DOUBLE);
} 

int
ncmpi_get_var_uchar_all(int ncid, int varid, unsigned char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_var_schar_all(int ncid, int varid, signed char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_var_text_all(int ncid, int varid, char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_var_short_all(int ncid, int varid, short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_var_int_all(int ncid, int varid, int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_INT);
} 

int
ncmpi_get_var_long_all(int ncid, int varid, long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_var_float_all(int ncid, int varid, float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_FLOAT);
} 

int
ncmpi_get_var_double_all(int ncid, int varid, double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nelems, MPI_DOUBLE);
} 

int
ncmpi_put_vara_uchar_all(int ncid, int varid,
                         const MPI_Offset start[], const MPI_Offset count[],
                         const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_vara_uchar(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_vara_schar_all(int ncid, int varid,
                         const MPI_Offset start[], const MPI_Offset count[],
                         const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_vara_schar(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_vara_text_all(int ncid, int varid,
                         const MPI_Offset start[], const MPI_Offset count[],
                         const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_vara_text(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_vara_short_all(int ncid, int varid,
                         const MPI_Offset start[], const MPI_Offset count[],
                         const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_vara_short(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_SHORT);
} 

int
ncmpi_put_vara_int_all(int ncid, int varid, 
                       const MPI_Offset start[], const MPI_Offset count[], 
		       const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;                                                                  

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count, 
			    (const void *)op, nelems, MPI_INT);
}

int
ncmpi_put_vara_int(int ncid, int varid, 
		const MPI_Offset start[], const MPI_Offset count[], 
		const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_INT);
}

int
ncmpi_put_vara_long_all(int ncid, int varid,
                       const MPI_Offset start[], const MPI_Offset count[],
                       const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_vara_long(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_vara_float_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count, 
			    (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_vara_float(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_vara_double_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vara_all(ncid, varid, start, count, 
                            (const void *)op, nelems, MPI_DOUBLE);
}

int
ncmpi_put_vara_double(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_DOUBLE);
}

int
ncmpi_get_vara_uchar_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    unsigned char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_vara_uchar(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    unsigned char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_vara_schar_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    signed char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_vara_schar(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    signed char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_vara_text_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_vara_text(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_vara_short_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    short *ip) {
 
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_vara_short(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    short *ip) {
 
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_SHORT);
} 

int
ncmpi_get_vara_int_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    int *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  
  return ncmpi_get_vara_all(ncid, varid, start, count, 
                            (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_vara_int(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_vara_long_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    long *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_vara_long(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_vara_float_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    float *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  
  return ncmpi_get_vara_all(ncid, varid, start, count, 
			    (void *)ip, nelems, MPI_FLOAT); 
}

int
ncmpi_get_vara_float(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_FLOAT);
}

int
ncmpi_get_vara_double_all(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vara_all(ncid, varid, start, count, 
			    (void *)ip, nelems, MPI_DOUBLE); 
}

int
ncmpi_get_vara_double(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_DOUBLE);
}

int
ncmpi_put_vars_uchar_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_vars_uchar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const unsigned char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_vars_schar_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_vars_schar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const signed char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_vars_text_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_vars_text(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_vars_short_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_vars_short(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_vars_int_all(int ncid, int varid,
                       const MPI_Offset start[],
                       const MPI_Offset count[],
                       const MPI_Offset stride[],
                       const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_INT);
}

int
ncmpi_put_vars_int(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
  
  nelems = 1; 
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
    
  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_INT);
} 

int
ncmpi_put_vars_long_all(int ncid, int varid,
                       const MPI_Offset start[],
                       const MPI_Offset count[],
                       const MPI_Offset stride[],
                       const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_vars_long(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   const long *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_vars_float_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_vars_float(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_vars_double_all(int ncid, int varid,
                          const MPI_Offset start[], 
                          const MPI_Offset count[],
                          const MPI_Offset stride[],
                          const double *op) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nelems, MPI_DOUBLE);

}

int
ncmpi_put_vars_double(int ncid, int varid,
                      const MPI_Offset start[],
                      const MPI_Offset count[],
                      const MPI_Offset stride[],
                      const double *op) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_DOUBLE);

}

int
ncmpi_get_vars_uchar_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         unsigned char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_vars_uchar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     unsigned char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_vars_schar_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         signed char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_vars_schar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     signed char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_vars_text_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_vars_text(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_vars_short_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_vars_short(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_vars_int_all(int ncid, int varid,
                       const MPI_Offset start[],
                       const MPI_Offset count[],
                       const MPI_Offset stride[],
                       int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_vars_int(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_vars_long_all(int ncid, int varid,
                       const MPI_Offset start[],
                       const MPI_Offset count[],
                       const MPI_Offset stride[],
                       long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_vars_long(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   long *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_vars_float_all(int ncid, int varid,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         const MPI_Offset stride[],
                         float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nelems, MPI_FLOAT);
}

int
ncmpi_get_vars_float(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_FLOAT);
}

int
ncmpi_get_vars_double_all(int ncid, int varid,
                          const MPI_Offset start[], 
			  const MPI_Offset count[],
			  const MPI_Offset stride[],
                          double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
			    (void *)ip, nelems, MPI_DOUBLE); 
}

int
ncmpi_get_vars_double(int ncid, int varid,
                      const MPI_Offset start[],
                      const MPI_Offset count[],
                      const MPI_Offset stride[],
                      double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_DOUBLE);
}

int
ncmpi_put_varm_uchar_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 const unsigned char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_varm_uchar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const unsigned char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_put_varm_schar_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 const signed char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_varm_schar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const signed char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_BYTE);
}

int
ncmpi_put_varm_text_all(int ncid, int varid,
			const MPI_Offset start[],
			const MPI_Offset count[],
			const MPI_Offset stride[],
			const MPI_Offset imap[],
			const char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_varm_text(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    const char *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_CHAR);
}

int
ncmpi_put_varm_short_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 const short *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_varm_short(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const short *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_SHORT);
}

int
ncmpi_put_varm_int_all(int ncid, int varid,
		       const MPI_Offset start[],
		       const MPI_Offset count[],
		       const MPI_Offset stride[],
		       const MPI_Offset imap[],
		       const int *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_INT);
}

int
ncmpi_put_varm_int(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   const int *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_INT);
}

int
ncmpi_put_varm_long_all(int ncid, int varid,
			const MPI_Offset start[],
			const MPI_Offset count[],
			const MPI_Offset stride[],
			const MPI_Offset imap[],
			const long *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_varm_long(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    const long *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_LONG);
}

int
ncmpi_put_varm_float_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 const float *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_varm_float(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const float *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_FLOAT);
}

int
ncmpi_put_varm_double_all(int ncid, int varid,
			  const MPI_Offset start[],
			  const MPI_Offset count[],
			  const MPI_Offset stride[],
			  const MPI_Offset imap[],
			  const double *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm_all(ncid, varid, start, count, stride, imap,
			    (const void *)op, nelems, MPI_DOUBLE);
}

int
ncmpi_put_varm_double(int ncid, int varid,
		      const MPI_Offset start[],
		      const MPI_Offset count[],
		      const MPI_Offset stride[],
		      const MPI_Offset imap[],
		      const double *op)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_put_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_DOUBLE);
}

int
ncmpi_get_varm_uchar_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 unsigned char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_varm_uchar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     unsigned char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR);
}

int
ncmpi_get_varm_schar_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 signed char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_varm_schar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     signed char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_BYTE);
}

int
ncmpi_get_varm_text_all(int ncid, int varid,
			const MPI_Offset start[],
			const MPI_Offset count[],
			const MPI_Offset stride[],
			const MPI_Offset imap[],
			char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_varm_text(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    char *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_CHAR);
}

int
ncmpi_get_varm_short_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 short *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_varm_short(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     short *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_SHORT);
}

int
ncmpi_get_varm_int_all(int ncid, int varid,
		       const MPI_Offset start[],
		       const MPI_Offset count[],
		       const MPI_Offset stride[],
		       const MPI_Offset imap[],
		       int *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_varm_int(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   int *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_INT);
}

int
ncmpi_get_varm_long_all(int ncid, int varid,
			const MPI_Offset start[],
			const MPI_Offset count[],
			const MPI_Offset stride[],
			const MPI_Offset imap[],
			long *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_varm_long(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    long *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_LONG);
}

int
ncmpi_get_varm_float_all(int ncid, int varid,
			 const MPI_Offset start[],
			 const MPI_Offset count[],
			 const MPI_Offset stride[],
			 const MPI_Offset imap[],
			 float *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_FLOAT);
}

int
ncmpi_get_varm_float(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     float *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_FLOAT);
}

int
ncmpi_get_varm_double_all(int ncid, int varid,
			  const MPI_Offset start[],
			  const MPI_Offset count[],
			  const MPI_Offset stride[],
			  const MPI_Offset imap[],
			  double *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm_all(ncid, varid, start, count, stride, imap,
			    (void *)ip, nelems, MPI_DOUBLE);
}

int
ncmpi_get_varm_double(int ncid, int varid,
		      const MPI_Offset start[],
		      const MPI_Offset count[],
		      const MPI_Offset stride[],
		      const MPI_Offset imap[],
		      double *ip)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_get_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_DOUBLE);
}

/* End {put,get}_var */

#ifdef ENABLE_NONBLOCKING
/* #################################################################### */
/* Begin non-blocking data access functions */

NCMPI_Request ncmpii_new_NCMPI_Request(int ncid, int varid, int ndims,
		const MPI_Offset *start, const MPI_Offset *count,
		const void *buf, MPI_Offset bufcount, MPI_Datatype datatype)
{
	NCMPI_Request request;
	int dim;

	request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));
	request->indep = 0;
	request->ncid = ncid;
	request->varid = varid;
	request->ndim = ndims;
	request->start = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	request->count = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	request->buf = (void *)buf;
	request->bufcount = bufcount;
	request->mpi_varatype = datatype;
	request->next_req = NULL;
	for (dim = 0; dim < ndims; dim++){
		request->start[dim]=start[dim];
		request->count[dim]=count[dim];
	}
	request->rw_flag = 1; 
	return request;
}

int ncmpii_free_NCMPI_Request(NCMPI_Request *request) 
{
	free ((*request)->start);
	free ((*request)->count);
	free (*request);
	return 0;
}

static void ncmpii_postwrite(void *xbuf, void *cbuf, void *buf) {
  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
}

static int ncmpii_postread(nc_type vartype,
		    void *xbuf, 
		    void *cbuf, 
		    MPI_Offset nelems, 
		    int cnelems, 
		    int iscontig_of_ptypes,
		    void *buf, 
		    int bufcount, 
		    MPI_Datatype datatype, 
		    MPI_Datatype ptype) 
{

  int status = NC_NOERR;

  if ( need_convert(vartype, ptype) ) {
    switch( vartype ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_getn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }
  } else if ( need_swap(vartype, ptype) ) {
if (cbuf == xbuf)
    in_swapn(xbuf, nelems, ncmpix_len_nctype(vartype));
else
    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(vartype));
  }

  if (!iscontig_of_ptypes) {
    status = ncmpii_data_repack(cbuf, cnelems, ptype,
                                (void *)buf, bufcount, datatype);
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  return status;
}

static void ncmpii_postmwrite(void *cbuf, void *lbuf, int iscontig_of_ptypes) {
  if (!iscontig_of_ptypes && lbuf != NULL)
    free(lbuf);
  if (cbuf != NULL)
    free(cbuf);
}

static int ncmpii_postmread(void *cbuf,
		     void *lbuf,
		     int cnelems,
		     int lnelems,
		     int iscontig_of_ptypes,
		     void *buf,
		     int bufcount,
                     MPI_Datatype datatype,
		     MPI_Datatype ptype,
                     MPI_Datatype imaptype)
{
  int status; 

  status = ncmpii_data_repack(cbuf, cnelems, ptype,
                              lbuf, 1, imaptype);
  if (status != NC_NOERR)
    return status;

  MPI_Type_free(&imaptype);

  if (!iscontig_of_ptypes) {
    status = ncmpii_data_repack(lbuf, lnelems, ptype,
                                (void *)buf, bufcount, datatype);
    if (lbuf != NULL)
      free(lbuf);
  }

  if (cbuf != NULL)
    free(cbuf);

  return status;
}

static int 
ncmpii_postprocess(NCMPI_Request *request) {
  int status = NC_NOERR;

  switch ((*request)->reqtype) {
    case NCMPI_REQTYPE_READ:
      status = ncmpii_postread((*request)->vartype,
			       (*request)->xbuf,
			       (*request)->cbuf,
			       (*request)->nelems,
			       (*request)->cnelems,
			       (*request)->iscontig_of_ptypes,
			       (*request)->buf,
			       (*request)->bufcount,
			       (*request)->datatype,
			       (*request)->ptype);
      MPI_Type_free(&((*request)->datatype));
      break;
    case NCMPI_REQTYPE_WRITE:
      ncmpii_postwrite((*request)->xbuf, (*request)->cbuf, (*request)->buf);
      break;
    case NCMPI_REQTYPE_MREAD:
      ncmpii_postmread((*request)->cbuf,
		       (*request)->lbuf,
		       (*request)->cnelems,
		       (*request)->lnelems,
		       (*request)->iscontig_of_ptypes,
		       (*request)->buf,
		       (*request)->bufcount,
		       (*request)->datatype,
		       (*request)->ptype,
		       (*request)->imaptype);
      MPI_Type_free(&((*request)->datatype));
      break;
    case NCMPI_REQTYPE_MWRITE:
      ncmpii_postmwrite((*request)->cbuf, 
			(*request)->lbuf, 
			(*request)->iscontig_of_ptypes);
      break;
    default:
      break;
  }
  if (status != NC_NOERR)
    return status;

  if ( (*request)->next_req != NCMPI_REQUEST_NULL ) {
    status = ncmpii_postprocess( &((*request)->next_req) );
    if (status != NC_NOERR)
      return status;
  }

  free(*request);
  *request = NCMPI_REQUEST_NULL;

  return status;
}
      

int
ncmpi_wait_one(NCMPI_Request *request) {
    int mpireturn = MPI_SUCCESS;

    if ((*request)->indep==1) {
        if (*request != NCMPI_REQUEST_NULL) {
            mpireturn = MPI_Wait(&((*request)->mpi_req), MPI_STATUS_IGNORE);
            ncmpii_postprocess(request);
        }
        if (mpireturn != MPI_SUCCESS)
            return NC_EFILE;
        else
            return NC_NOERR;
    } else {
        return (ncmpi_coll_wait(*request));
    }
}

int
ncmpi_waitall(int count, NCMPI_Request array_of_requests[]) {
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  if(array_of_requests[0]->indep==1) {
  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Waitall(count, array_of_mpireqs, MPI_STATUSES_IGNORE);
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      ncmpii_postprocess(array_of_requests+i);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
    }
  } else {
    return (ncmpi_coll_waitall(count, array_of_requests));
  }
}

int
ncmpi_wait_all(int count, NCMPI_Request array_of_requests[], int array_of_statuses[]) {
  int i,j;
  int ncid;
  int nvars;
  int *varids;
  MPI_Offset **starts;
  MPI_Offset **counts;
  void **buf;
  MPI_Offset *bufcount;
  MPI_Datatype *datatype;
  int ndim;
  int status;

  ncid = array_of_requests[0]->ncid;
  nvars = count;
  varids = (int *)malloc(count*sizeof(int));
  starts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  counts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  buf   = (void **)malloc(count*sizeof(void *));
  bufcount = (MPI_Offset *)malloc(count*sizeof(MPI_Offset));
  datatype = (MPI_Datatype *)malloc(count*sizeof(MPI_Datatype));

  for (i=0; i<count; i++){
        ndim = array_of_requests[i]->ndim;
        varids[i]=array_of_requests[i]->varid;
        starts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));
        counts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));
        for (j=0; j<ndim; j++){
           starts[i][j]=array_of_requests[i]->start[j];
           counts[i][j]=array_of_requests[i]->count[j];
        }
        buf[i] = array_of_requests[i]->buf;
        bufcount[i] = array_of_requests[i]->bufcount;
        datatype[i] = array_of_requests[i]->mpi_varatype;
  }
  if (array_of_requests[0]->rw_flag == 1)
  status = ncmpi_put_mvara_all(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);

  if (array_of_requests[0]->rw_flag == 0)
  status = ncmpi_get_mvara_all(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);

  for (i=0; i<count; i++){
	array_of_statuses[i] = status;
	ncmpii_free_NCMPI_Request(&(array_of_requests[i]));

  }
  free(starts);
  free(varids);
  free(counts);
  free(buf);
  free(bufcount);
  free(datatype);
  return NC_NOERR;
}


int
ncmpi_wait(int count, NCMPI_Request array_of_requests[], int array_of_statuses[]) {
  int i,j;
  int ncid;
  int nvars;
  int *varids;
  MPI_Offset **starts;
  MPI_Offset **counts;
  void **buf;
  MPI_Offset *bufcount;
  MPI_Datatype *datatype;
  int ndim;
  int status;

  ncid = array_of_requests[0]->ncid;
  nvars = count;
  varids = (int *)malloc(count*sizeof(int));
  starts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  counts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  buf   = (void **)malloc(count*sizeof(void *));
  bufcount = (MPI_Offset *)malloc(count*sizeof(MPI_Offset));
  datatype = (MPI_Datatype *)malloc(count*sizeof(MPI_Datatype));
 
  status = ncmpi_begin_indep_data(ncid);

  for (i=0; i<count; i++){
        ndim = array_of_requests[i]->ndim;
        varids[i]=array_of_requests[i]->varid;
        starts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));
        counts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));
        for (j=0; j<ndim; j++){
           starts[i][j]=array_of_requests[i]->start[j];
           counts[i][j]=array_of_requests[i]->count[j];
        }
        buf[i] = array_of_requests[i]->buf;
        bufcount[i] = array_of_requests[i]->bufcount;
        datatype[i] = array_of_requests[i]->mpi_varatype;
  }
  if (array_of_requests[0]->rw_flag == 1)
  status = ncmpi_put_mvara(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);

  if (array_of_requests[0]->rw_flag == 0)
  status = ncmpi_get_mvara(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);

  for (i=0; i<count; i++){
        array_of_statuses[i] = status;
  }
  status = ncmpi_end_indep_data(ncid);
  free(starts);
  free(varids);
  free(counts);
  free(buf);
  free(bufcount);
  free(datatype);
  return NC_NOERR;
}


int
ncmpi_waitany(int count, NCMPI_Request array_of_requests[], int *index) {
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Waitany(count, array_of_mpireqs, index, MPI_STATUS_IGNORE);
  if (array_of_requests[*index] != NCMPI_REQUEST_NULL)
    ncmpii_postprocess(array_of_requests + *index);

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_waitsome(int count, NCMPI_Request array_of_requests[],
	       int *outcount, int array_of_indices[]) 
{
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Waitsome(count, array_of_mpireqs, 
			   outcount, array_of_indices, MPI_STATUSES_IGNORE);
  for (i=0; i<*outcount; i++) {
    if (array_of_requests[array_of_indices[i]] != NCMPI_REQUEST_NULL)
      ncmpii_postprocess(array_of_requests + array_of_indices[i]);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_test(NCMPI_Request *request, int *flag) {
  int mpireturn = MPI_SUCCESS;
  if (*request != NCMPI_REQUEST_NULL) {
    mpireturn = MPI_Test(&((*request)->mpi_req), flag, MPI_STATUS_IGNORE);
    if (*flag) 
      ncmpii_postprocess(request);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_testall(int count, NCMPI_Request array_of_requests[], int *flag) {
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL) 
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Testall(count, array_of_mpireqs, flag, MPI_STATUSES_IGNORE);
  if (*flag) {
    for (i=0; i<count; i++) {
      if (array_of_requests[i] != NCMPI_REQUEST_NULL) 
	ncmpii_postprocess(array_of_requests+i);
    }
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_testany(int count, NCMPI_Request array_of_requests[], 
	      int *index, int *flag)
{
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Testany(count, array_of_mpireqs, 
			  index, flag, MPI_STATUS_IGNORE);
  if (*flag && array_of_requests[*index] != NCMPI_REQUEST_NULL)
    ncmpii_postprocess(array_of_requests + *index);

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_testsome(int count, NCMPI_Request array_of_requests[],
	       int *outcount, int array_of_indices[])  
{
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (MPI_Request *)malloc(count * sizeof(int));
  for (i=0; i<count; i++) {
    if (array_of_requests[i] != NCMPI_REQUEST_NULL)
      array_of_mpireqs[i] = array_of_requests[i]->mpi_req;
    else
      array_of_mpireqs[i] = MPI_REQUEST_NULL;
  }
  mpireturn = MPI_Testsome(count, array_of_mpireqs,
                           outcount, array_of_indices, MPI_STATUSES_IGNORE);
  for (i=0; i<*outcount; i++) {
    if (array_of_requests[array_of_indices[i]] != NCMPI_REQUEST_NULL)
      ncmpii_postprocess(array_of_requests + array_of_indices[i]);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_request_get_status(NCMPI_Request request, int *flag) {
  int mpireturn = MPI_SUCCESS;
  if (request != NCMPI_REQUEST_NULL)
    mpireturn = MPI_Request_get_status(request->mpi_req, flag, 
				       MPI_STATUS_IGNORE);
  else
    *flag = 1;

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_request_free(NCMPI_Request *request) {
  int mpireturn = MPI_SUCCESS;
  if (*request != NCMPI_REQUEST_NULL) {
    mpireturn = MPI_Request_free( &((*request)->mpi_req) );
    ncmpii_postprocess(request);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int
ncmpi_cancel(NCMPI_Request *request) {
  int mpireturn = MPI_SUCCESS;
  if (*request != NCMPI_REQUEST_NULL) 
    mpireturn = MPI_Cancel(&((*request)->mpi_req));

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}


int
ncmpi_iput_var1(int ncid, int varid,
               const MPI_Offset index[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes, offset;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status; 

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  if (nelems != cnelems) {
    if (warning == NC_NOERR) 
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems*varp->xsz; /* account for file bytes */
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* accessing one variable element need not set the file view.
     just find the file offset and use MPI-IO call with explicit offset */

  status = get_offset(ncp, varp, index, &offset);
  if (status != NC_NOERR) return status;

  if (!iscontig_of_ptypes) {
  
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype, 
				cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
   
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */
 
  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */
  
    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, 1, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, 1, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, 1, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, 1, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, 1, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, 1, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite_at(ncp->nciop->independent_fh, offset, xbuf, nbytes,
			         MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iwrite_at");
        status = NC_EWRITE;
  }
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  if ( status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = index[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_iget_var1(int ncid, int varid,
               const MPI_Offset index[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request) 
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes, offset;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status; 

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems*varp->xsz; /* account for file bytes */
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* accessing one variable element need not set the file view.
     just find the file offset and use MPI-IO call with explicit offset */

  status = get_offset(ncp, varp, index, &offset);
  if (status != NC_NOERR) return status;

  if (!iscontig_of_ptypes) {
  
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) || need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iread_at(ncp->nciop->independent_fh, offset, xbuf, 
			        nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_read_at");
        status = NC_EREAD;
  }
 
  (*request)->reqtype = NCMPI_REQTYPE_READ;
  (*request)->vartype = varp->type;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->nelems = nelems;
  (*request)->cnelems = cnelems;
  (*request)->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->buf = (void *)buf;
  (*request)->bufcount = bufcount;
  MPI_Type_dup(datatype, &((*request)->datatype));
  (*request)->ptype = ptype;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_iput_var(int ncid, int varid, 
	       const void *buf, MPI_Offset bufcount, MPI_Datatype datatype,
	       NCMPI_Request *request) 
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  if (varp->ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (varp->ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems * varp->xsz;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* set the mpi file view */
 
  status = set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
  if(status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iwrite");
        status = NC_EWRITE;
  }
 
  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    if (varp->ndims > 1)
      newnumrecs = nelems / varp->dsizes[1];
    else
      newnumrecs = nelems;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_iget_var(int ncid, int varid, 
	       void *buf, MPI_Offset bufcount, MPI_Datatype datatype, 
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  if (varp->ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (varp->ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = nelems * varp->xsz;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ||  need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iread(ncp->nciop->independent_fh, xbuf, 
			     nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iread");
        status = NC_EREAD;
  }
 
  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_READ;
  (*request)->vartype = varp->type;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->nelems = nelems;
  (*request)->cnelems = cnelems;
  (*request)->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->buf = (void *)buf;
  (*request)->bufcount = bufcount;
  MPI_Type_dup(datatype, &((*request)->datatype));
  (*request)->ptype = ptype;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_iput_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL,  *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
  if(NC_indep(ncp)) {
    status = ncmpii_NC_check_id(ncid, &ncp);
    if(status != NC_NOERR)
      return status;
 
    if(NC_readonly(ncp))
      return NC_EPERM;
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp, start, count, 0);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, ptype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iwrite");
        return NC_EWRITE;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;
    (*request)->indep = 1;

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = start[0] + count[0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
  } else {
  return ncmpi_iput_vara_all(ncid, varid, start, count,
                        buf, bufcount, datatype, request);
  }
} 

int
ncmpi_iget_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;

  /* set the mpi file view */
 
  status = set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp, start, count, 1);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }
 
  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ||  need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iread(ncp->nciop->independent_fh, xbuf, 
			     nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iread");
        status = NC_EREAD;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_READ;
  (*request)->indep = 1;
  (*request)->vartype = varp->type;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->nelems = nelems;
  (*request)->cnelems = cnelems;
  (*request)->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->buf = (void *)buf;
  (*request)->bufcount = bufcount;
  MPI_Type_dup(datatype, &((*request)->datatype));
  (*request)->ptype = ptype;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_iput_vars(int ncid, int varid,
               const MPI_Offset start[], 
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
			        varp, start, count, stride, 0);
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                cbuf, cnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ) {

    /* allocate new buffer */

    xbuf = (void *)malloc(nbytes);

    /* automatic numeric datatype conversion */

    switch( varp->type ) {
      case NC_BYTE:
         status = x_putn_schar(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_SHORT:
         status = x_putn_short(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_INT:
         status = x_putn_int(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_FLOAT:
         status = x_putn_float(xbuf, cbuf, cnelems, ptype);
         break;
      case NC_DOUBLE:
         status = x_putn_double(xbuf, cbuf, cnelems, ptype);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp->type, datatype) ) {

    in_swapn(cbuf, nelems, ncmpix_len_nctype(varp->type));
    xbuf = (void *)cbuf;

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iwrite");
        status = NC_EWRITE;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    MPI_Offset newnumrecs;
    newnumrecs = start[0] + (count[0] - 1) * stride[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_iget_vars(int ncid, int varid,
               const MPI_Offset start[], 
	       const MPI_Offset count[],
               const MPI_Offset stride[],
               void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset nelems, cnelems, nbytes;
  int el_size;
  int mpireturn;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &cnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp->type, ptype) )
    return NC_ECHAR;

  cnelems *= bufcount;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    nelems *= count[dim];
  }

  if (nelems != cnelems) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);
  }

  nbytes = varp->xsz * nelems;
  if (nbytes < 0)
    return NC_ENEGATIVECNT;
 
  /* set the mpi file view */
 
  status = set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
				varp, start, count, stride, 1); 
  if(status != NC_NOERR)
    return status;
 
  if (!iscontig_of_ptypes) {
 
    /* account for derived datatype: allocate the contiguous buffer */
 
    cbuf = (void *)malloc( cnelems * el_size );
 
  } else {
 
    cbuf = (void *)buf;
 
  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp->type, ptype) ||  need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iread(ncp->nciop->independent_fh, xbuf, 
			     nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_iread");
        status = NC_EREAD;
  }
 
  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  (*request)->reqtype = NCMPI_REQTYPE_READ;
  (*request)->vartype = varp->type;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->nelems = nelems;
  (*request)->cnelems = cnelems;
  (*request)->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->buf = (void *)buf;
  (*request)->bufcount = bufcount;
  MPI_Type_dup(datatype, &((*request)->datatype));
  (*request)->ptype = ptype;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_iput_varm(int ncid, int varid,
	       const MPI_Offset start[],
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
	       const MPI_Offset imap[],
	       const void *buf, MPI_Offset bufcount,
	       MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems; 
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_iput_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype, request);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_iput_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype, request);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_iput_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype, request);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {

    /* handling for derived datatype: pack into a contiguous buffer */

    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      return status;

  } else {

    lbuf = (void *)buf;

  }

  if (count[dim] < 0)
    return NC_ENEGATIVECNT;
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
                  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen*count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
  
#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, imap[dim]*el_size,
                            imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  /* layout lbuf to cbuf based on imap */
  status = ncmpii_data_repack(lbuf, 1, imaptype,
			      cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    return status;

  MPI_Type_free(&imaptype);

  status = ncmpi_iput_vars(ncid, varid, start, count, stride,
                          cbuf, cnelems, ptype, request);

  (*request)->next_req = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));
  (*request)->next_req->reqtype = NCMPI_REQTYPE_MWRITE;
  (*request)->next_req->cbuf = cbuf;
  (*request)->next_req->lbuf = lbuf;
  (*request)->next_req->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->next_req->next_req = NCMPI_REQUEST_NULL;

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_iget_varm(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   void *buf, MPI_Offset bufcount,
		   MPI_Datatype datatype,
		   NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  MPI_Offset ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  MPI_Offset lnelems, cnelems; 
  int el_size;
  MPI_Datatype ptype, tmptype, imaptype;
  int isderived, iscontig_of_ptypes;
  int imap_contig_blocklen;

  if (imap == NULL) {
    /* no mapping, same as vars */
    return ncmpi_iget_vars(ncid, varid, start, count, stride,
                          buf, bufcount, datatype, request);
  }

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims == 0) {
    /* reduced to scalar var, only one value at one fixed place */
    return ncmpi_iget_vars(ncid, varid, start, count, stride,
			  buf, bufcount, datatype, request);
  }

  imap_contig_blocklen = 1;
  dim = ndims;
  /* test each dim's contiguity until the 1st non-contiguous dim is reached */
  while ( --dim>=0 && imap_contig_blocklen==imap[dim] ) {
    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
    imap_contig_blocklen *= count[dim];
  }

  if (dim == -1) {
    /* imap is a contiguous layout */
    return ncmpi_iget_vars(ncid, varid, start, count, stride,
			  buf, bufcount, datatype, request);
  } /* else imap gives non-contiguous layout, and need pack/unpack */

  status = ncmpii_dtype_decode(datatype, &ptype, &el_size,
			       &lnelems, &isderived, &iscontig_of_ptypes);
  if (status != NC_NOERR)
    return status;

  if (!iscontig_of_ptypes) {
 
    /* handling for derived datatype: pack into a contiguous buffer */
 
    lnelems *= bufcount;
    lbuf = (void *)malloc( lnelems*el_size );
    status = ncmpii_data_repack((void *)buf, bufcount, datatype,
                                lbuf, lnelems, ptype);
    if (status != NC_NOERR)
      return status;
 
  } else {
 
    lbuf = (void *)buf;
 
  }

  if (count[dim] < 0)
    return NC_ENEGATIVECNT;
  MPI_Type_vector(count[dim], imap_contig_blocklen, imap[dim],
		  ptype, &imaptype);
  MPI_Type_commit(&imaptype);
  cnelems = imap_contig_blocklen * count[dim];
  for (dim--; dim>=0; dim--) {

    if (count[dim] < 0)
      return NC_ENEGATIVECNT;
  
#if (MPI_VERSION < 2)
    MPI_Type_hvector(count[dim], 1, imap[dim]*el_size, imaptype, &tmptype);
#else
    MPI_Type_create_hvector(count[dim], 1, (MPI_Aint)imap[dim]*el_size, 
			    imaptype, &tmptype);
#endif
    MPI_Type_free(&imaptype);
    MPI_Type_commit(&tmptype);
    imaptype = tmptype;
    cnelems *= count[dim];

  }

  cbuf = (void *)malloc(cnelems*el_size);

  status = ncmpi_iget_vars(ncid, varid, start, count, stride, 
			  cbuf, cnelems, ptype, request);
  if (status != NC_NOERR) {
    if (status == NC_ERANGE && warning == NC_NOERR) 
      warning = status; /* to satisfy the nc_test logic */
    else
      return status;
  }

  (*request)->next_req = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));
  (*request)->next_req->reqtype = NCMPI_REQTYPE_MREAD;
  (*request)->next_req->cbuf = cbuf;
  (*request)->next_req->lbuf = lbuf;
  (*request)->next_req->cnelems = cnelems;
  (*request)->next_req->lnelems = lnelems;
  (*request)->next_req->iscontig_of_ptypes = iscontig_of_ptypes;
  (*request)->next_req->buf = (void *)buf;
  (*request)->next_req->bufcount = bufcount;
  MPI_Type_dup(datatype, &((*request)->next_req->datatype));
  (*request)->next_req->ptype = ptype;
  (*request)->next_req->imaptype = imaptype;
    
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_iput_var1_uchar(int ncid, int varid,
                     const MPI_Offset index[],
                     const unsigned char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iput_var1_schar(int ncid, int varid,
                     const MPI_Offset index[],
                     const signed char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_BYTE, request);
}

int
ncmpi_iput_var1_text(int ncid, int varid,
                     const MPI_Offset index[],
                     const char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_CHAR, request);
}


int
ncmpi_iput_var1_short(int ncid, int varid,
                     const MPI_Offset index[],
		     const short *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iput_var1(ncid, varid, index, 
                        (const void *)op, 1, MPI_SHORT, request); 
}

int
ncmpi_iput_var1_int(int ncid, int varid,
                   const MPI_Offset index[],
                   const int *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_INT, request);
}

int
ncmpi_iput_var1_long(int ncid, int varid,
                   const MPI_Offset index[],
                   const long *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_LONG, request);
}

int
ncmpi_iput_var1_float(int ncid, int varid,
                     const MPI_Offset index[],
                     const float *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_FLOAT, request); 
}
 
int
ncmpi_iput_var1_double(int ncid, int varid,
                      const MPI_Offset index[],
                      const double *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iput_var1(ncid, varid, index,
                        (const void *)op, 1, MPI_DOUBLE, request);
}

int
ncmpi_iget_var1_uchar(int ncid, int varid,
                     const MPI_Offset index[],
                     unsigned char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iget_var1_schar(int ncid, int varid,
                     const MPI_Offset index[],
                     signed char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_BYTE, request);
}

int
ncmpi_iget_var1_text(int ncid, int varid,
                     const MPI_Offset index[],
                     char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_CHAR, request);
}

int
ncmpi_iget_var1_short(int ncid, int varid,
                     const MPI_Offset index[],
                     short *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_SHORT, request); 
}
 
int
ncmpi_iget_var1_int(int ncid, int varid,
                   const MPI_Offset index[],
                   int *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_INT, request);  
}
 
int
ncmpi_iget_var1_long(int ncid, int varid,
                   const MPI_Offset index[],
                   long *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_LONG, request); 
}

int
ncmpi_iget_var1_float(int ncid, int varid,
                     const MPI_Offset index[],
                     float *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_FLOAT, request);  
}
 
int
ncmpi_iget_var1_double(int ncid, int varid,
                      const MPI_Offset index[],
                      double *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_iget_var1(ncid, varid, index,
                        (void *)ip, 1, MPI_DOUBLE, request);  
} 

int
ncmpi_iput_var_uchar(int ncid, int varid, const unsigned char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iput_var_schar(int ncid, int varid, const signed char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_BYTE, request);
}


int
ncmpi_iput_var_text(int ncid, int varid, const char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */
  
  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;
  
  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_CHAR, request);
}

int
ncmpi_iput_var_short(int ncid, int varid, const short *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_SHORT, request);
}

int
ncmpi_iput_var_int(int ncid, int varid, const int *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_INT, request);
} 

int
ncmpi_iput_var_long(int ncid, int varid, const long *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_LONG, request);
}

int
ncmpi_iput_var_float(int ncid, int varid, const float *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_FLOAT, request);
} 

int
ncmpi_iput_var_double(int ncid, int varid, const double *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iput_var(ncid, varid, (const void *)op, nelems, MPI_DOUBLE, request);
} 

int
ncmpi_iget_var_uchar(int ncid, int varid, unsigned char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iget_var_schar(int ncid, int varid, signed char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_BYTE, request);
}

int
ncmpi_iget_var_text(int ncid, int varid, char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_CHAR, request);
}

int
ncmpi_iget_var_short(int ncid, int varid, short *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_SHORT, request);
}

int
ncmpi_iget_var_int(int ncid, int varid, int *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_INT, request);
} 

int
ncmpi_iget_var_long(int ncid, int varid, long *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */

  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_LONG, request);
}

int
ncmpi_iget_var_float(int ncid, int varid, float *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_FLOAT, request);
} 

int
ncmpi_iget_var_double(int ncid, int varid, double *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  /* Begin modification 20030311 to fix bug of 0-dimensional variables */

  if (ndims == 0)
    nelems = 1;
  else if (!IS_RECVAR(varp))
    nelems = varp->dsizes[0];
  else if (ndims > 1)
    nelems = ncp->numrecs * varp->dsizes[1];
  else
    nelems = ncp->numrecs;

  /* End modification 20030311 */
 
  return ncmpi_iget_var(ncid, varid, (void *)ip, nelems, MPI_DOUBLE, request);
} 

int
ncmpi_iput_vara_uchar(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const unsigned char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iput_vara_schar(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const signed char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_BYTE, request);
}

int
ncmpi_iput_vara_text(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_CHAR, request);
}

int
ncmpi_iput_vara_short(int ncid, int varid,
                     const MPI_Offset start[], const MPI_Offset count[],
                     const short *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_SHORT, request);
} 

int
ncmpi_iput_vara_int(int ncid, int varid, 
		const MPI_Offset start[], const MPI_Offset count[], 
		const int *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_INT, request);
}

int
ncmpi_iput_vara_long(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const long *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_LONG, request);
}

int
ncmpi_iput_vara_float(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const float *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_FLOAT, request);
}

int
ncmpi_iput_vara_double(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                const double *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iput_vara(ncid, varid, start, count,
                        (const void *)op, nelems, MPI_DOUBLE, request);
}

int
ncmpi_iget_vara_uchar(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    unsigned char *ip,
		     NCMPI_Request *request)
{

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iget_vara_schar(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    signed char *ip,
		     NCMPI_Request *request)
{

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_BYTE, request);
}

int
ncmpi_iget_vara_text(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    char *ip,
		     NCMPI_Request *request)
{

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_CHAR, request);
}

int
ncmpi_iget_vara_short(int ncid, int varid,
                    const MPI_Offset start[], const MPI_Offset count[],
                    short *ip,
		     NCMPI_Request *request)
{
 
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_SHORT, request);
} 

int
ncmpi_iget_vara_int(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                int *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_INT, request);
}

int
ncmpi_iget_vara_long(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                long *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_LONG, request);
}

int
ncmpi_iget_vara_float(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                float *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_FLOAT, request);
}

int
ncmpi_iget_vara_double(int ncid, int varid,
                const MPI_Offset start[], const MPI_Offset count[],
                double *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
 
  return ncmpi_iget_vara(ncid, varid, start, count,
                        (void *)ip, nelems, MPI_DOUBLE, request);
}

int
ncmpi_iput_vars_uchar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const unsigned char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iput_vars_schar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const signed char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_BYTE, request);
}

int
ncmpi_iput_vars_text(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_CHAR, request);
}

int
ncmpi_iput_vars_short(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const short *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_SHORT, request);
}

int
ncmpi_iput_vars_int(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   const int *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
  
  nelems = 1; 
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
    
  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_INT, request);
} 

int
ncmpi_iput_vars_long(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   const long *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_LONG, request);
}

int
ncmpi_iput_vars_float(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     const float *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_FLOAT, request);
}

int
ncmpi_iput_vars_double(int ncid, int varid,
                      const MPI_Offset start[],
                      const MPI_Offset count[],
                      const MPI_Offset stride[],
                      const double *op,
		     NCMPI_Request *request)
{

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_vars(ncid, varid, start, count, stride,
                        (const void *)op, nelems, MPI_DOUBLE, request);

}

int
ncmpi_iget_vars_uchar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     unsigned char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iget_vars_schar(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     signed char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_BYTE, request);
}

int
ncmpi_iget_vars_text(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_CHAR, request);
}

int
ncmpi_iget_vars_short(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     short *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_SHORT, request);
}

int
ncmpi_iget_vars_int(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   int *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_INT, request);
}

int
ncmpi_iget_vars_long(int ncid, int varid,
                   const MPI_Offset start[],
                   const MPI_Offset count[],
                   const MPI_Offset stride[],
                   long *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_LONG, request);
}

int
ncmpi_iget_vars_float(int ncid, int varid,
                     const MPI_Offset start[],
                     const MPI_Offset count[],
                     const MPI_Offset stride[],
                     float *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_FLOAT, request);
}

int
ncmpi_iget_vars_double(int ncid, int varid,
                      const MPI_Offset start[],
                      const MPI_Offset count[],
                      const MPI_Offset stride[],
                      double *ip,
		     NCMPI_Request *request)
{

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_vars(ncid, varid, start, count, stride,
                        (void *)ip, nelems, MPI_DOUBLE, request);
}

int
ncmpi_iput_varm_uchar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const unsigned char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iput_varm_schar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const signed char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_BYTE, request);
}

int
ncmpi_iput_varm_text(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    const char *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_CHAR, request);
}

int
ncmpi_iput_varm_short(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const short *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_SHORT, request);
}

int
ncmpi_iput_varm_int(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   const int *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_INT, request);
}

int
ncmpi_iput_varm_long(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    const long *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  long int nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_LONG, request);
}

int
ncmpi_iput_varm_float(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     const float *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_FLOAT, request);
}

int
ncmpi_iput_varm_double(int ncid, int varid,
		      const MPI_Offset start[],
		      const MPI_Offset count[],
		      const MPI_Offset stride[],
		      const MPI_Offset imap[],
		      const double *op,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iput_varm(ncid, varid, start, count, stride, imap,
                        (const void *)op, nelems, MPI_DOUBLE, request);
}

int
ncmpi_iget_varm_uchar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     unsigned char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_UNSIGNED_CHAR, request);
}

int
ncmpi_iget_varm_schar(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     signed char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_BYTE, request);
}

int
ncmpi_iget_varm_text(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    char *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_CHAR, request);
}

int
ncmpi_iget_varm_short(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     short *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_SHORT, request);
}

int
ncmpi_iget_varm_int(int ncid, int varid,
		   const MPI_Offset start[],
		   const MPI_Offset count[],
		   const MPI_Offset stride[],
		   const MPI_Offset imap[],
		   int *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_INT, request);
}

int
ncmpi_iget_varm_long(int ncid, int varid,
		    const MPI_Offset start[],
		    const MPI_Offset count[],
		    const MPI_Offset stride[],
		    const MPI_Offset imap[],
		    long *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_LONG, request);
}

int
ncmpi_iget_varm_float(int ncid, int varid,
		     const MPI_Offset start[],
		     const MPI_Offset count[],
		     const MPI_Offset stride[],
		     const MPI_Offset imap[],
		     float *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_FLOAT, request);
}

int
ncmpi_iget_varm_double(int ncid, int varid,
		      const MPI_Offset start[],
		      const MPI_Offset count[],
		      const MPI_Offset stride[],
		      const MPI_Offset imap[],
		      double *ip,
		     NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  MPI_Offset nelems;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];

  return ncmpi_iget_varm(ncid, varid, start, count, stride, imap,
                        (void *)ip, nelems, MPI_DOUBLE, request);
}

static int
set_vara_fileview_all(NC* ncp, MPI_File *mpifh, int nvars, NC_var **varp, MPI_Offset *start[], MPI_Offset *count[], int getnotput) {

  MPI_Offset *offset;
  int status;
  int dim, *ndims;
  int *shape = NULL, *subcount = NULL, *substart = NULL; /* all in bytes */
  MPI_Datatype rectype;
  MPI_Datatype *filetype;
  MPI_Datatype full_filetype;
  int mpireturn;
  int i;
  int *blocklen;

  offset = (MPI_Offset *)malloc(sizeof(MPI_Offset)*nvars);
  ndims = (int *)malloc(sizeof(int)*nvars);
  blocklen = (int *)malloc(sizeof(int)*nvars);
  filetype = (MPI_Datatype *)malloc(sizeof(MPI_Datatype)*nvars);  

  for (i=0; i<nvars; i++){
  offset[i] = varp[i]->begin;
  ndims[i] = varp[i]->ndims;

  /* New coordinate/edge check to fix NC_EINVALCOORDS bug */
  status = NCedgeck(ncp, varp[i], start[i], count[i]);
  if( (status != NC_NOERR) ||
      (getnotput && IS_RECVAR(varp[i]) && *start[i] + *count[i] > NC_get_numrecs(ncp)) ) 
  {
    status = NCcoordck(ncp, varp[i], start[i]);
    if (status != NC_NOERR)
      return status;
    else
      return NC_EEDGE;
  }

/* Removed to fix NC_EINVALCOORDS bug

  status = NCcoordck(ncp, varp, start);
  if (status != NC_NOERR)
    return status;

  status = NCedgeck(ncp, varp, start, count);
  if(status != NC_NOERR)
    return status;

  if (getnotput && IS_RECVAR(varp) && *start + *count > NC_get_numrecs(ncp))
    return NC_EEDGE;
*/

  if (ndims[i] == 0) {

    /* scalar variable */
    filetype[i] = MPI_BYTE;

  } else {

    /* if ndims == 0, all below pointers would be null */

    shape = (int *)malloc(sizeof(int) * ndims[i]);
    subcount = (int *)malloc(sizeof(int) * ndims[i]);
    substart = (int *)malloc(sizeof(int) * ndims[i]);

    dim = 0;
    while (dim < ndims[i] && count[i][dim] > 0) dim++;

    if (dim < ndims[i]) {

      /* 0 size data */
      filetype[i] = MPI_BYTE;

    } else {
      if (IS_RECVAR(varp[i])){ 
        subcount[0] = 1;
//        subcount[0] = count[i][0];
//        substart[0] = start[i][0];
        substart[0] = 0;
        shape[0] = subcount[0];

        if (ncp->recsize <= varp[i]->len) {
    
          /* the only record variable */

          if (varp[i]->ndims == 1) {
            shape[0] *= varp[i]->xsz;
	    subcount[0] *= varp[i]->xsz;
          } else {
	    for (dim = 1; dim < ndims[i]-1; dim++) {
              shape[dim] = varp[i]->shape[dim];
              subcount[dim] = count[i][dim];
              substart[dim] = start[i][dim];
	    }
	    shape[dim] = varp[i]->xsz * varp[i]->shape[dim];
	    subcount[dim] = varp[i]->xsz * count[i][dim];
	    substart[dim] = varp[i]->xsz * start[i][dim];
          }
          offset[i] += start[i][0] * ncp->recsize;

          MPI_Type_create_subarray(ndims[i], shape, subcount, substart, 
				    MPI_ORDER_C, MPI_BYTE, &filetype[i]); 
    
          MPI_Type_commit(&filetype[i]);
        } else {
          /* more than one record variables */

          offset[i] += start[i][0] * ncp->recsize;
          if (varp[i]->ndims == 1) {
#if (MPI_VERSION < 2)
	    MPI_Type_hvector(subcount[0], varp[i]->xsz, ncp->recsize,
			    MPI_BYTE, &filetype[i]);
#else
	    MPI_Type_create_hvector(subcount[0], varp[i]->xsz, ncp->recsize,
				    MPI_BYTE, &filetype[i]);
#endif
	    MPI_Type_commit(&filetype[i]);

          } else {
            for (dim = 1; dim < ndims[i]-1; dim++) {
              shape[dim] = varp[i]->shape[dim];
              subcount[dim] = count[i][dim];
              substart[dim] = start[i][dim];
            }
            shape[dim] = varp[i]->xsz * varp[i]->shape[dim];
            subcount[dim] = varp[i]->xsz * count[i][dim];
            substart[dim] = varp[i]->xsz * start[i][dim];

	    MPI_Type_create_subarray(ndims[i]-1, shape+1, subcount+1, substart+1,
				     MPI_ORDER_C, MPI_BYTE, &rectype);
	    MPI_Type_commit(&rectype);
#if (MPI_VERSION < 2)
	    MPI_Type_hvector(subcount[0], 1, ncp->recsize, rectype, &filetype[i]);
#else
	    MPI_Type_create_hvector(subcount[0], 1, ncp->recsize, rectype, &filetype[i]);
#endif
	    MPI_Type_commit(&filetype[i]);
	    MPI_Type_free(&rectype);
          }
        }
      } else {

        /* non record variable */

        for (dim = 0; dim < ndims[i]-1; dim++ ) {
          shape[dim] = varp[i]->shape[dim];
          subcount[dim] = count[i][dim];
          substart[dim] = start[i][dim];
        }

        shape[dim] = varp[i]->xsz * varp[i]->shape[dim];
        subcount[dim] = varp[i]->xsz * count[i][dim];
        substart[dim] = varp[i]->xsz * start[i][dim];

        MPI_Type_create_subarray(ndims[i], shape, subcount, substart, 
		         MPI_ORDER_C, MPI_BYTE, &filetype[i]); 

        MPI_Type_commit(&filetype[i]);
      }
    }
  }
   if (IS_RECVAR(varp[i])&&(ncp->recsize <= varp[i]->len)) {
        blocklen[i]=count[i][0];
    } else {
        blocklen[i]=1;
    }

//    blocklen[i]=1;
    free(shape);
    free(subcount);
    free(substart);
  }/*end for loop*/

  MPI_Type_create_struct(nvars, blocklen, (MPI_Aint *)offset, filetype, &full_filetype);
  MPI_Type_commit(&full_filetype);
  
 // mpireturn = MPI_File_set_view(*mpifh, offset[0], MPI_BYTE, 
  mpireturn = MPI_File_set_view(*mpifh, 0, MPI_BYTE, 
		    full_filetype, "native", MPI_INFO_NULL);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_set_view");
        return NC_EFILE;
  }

  for (i=0; i<nvars; i++){
  if (ndims[i] > 0) {
    if (filetype[i] != MPI_BYTE) {
      MPI_Type_free(&filetype[i]);
    }	
  }
  }/*end for loop*/
  if (full_filetype != MPI_BYTE)
      MPI_Type_free(&full_filetype);
  free(offset);
  free(ndims);
  free(blocklen);
  free(filetype);
  return NC_NOERR;
}


static int
ncmpi_put_mvara_all_nonrecord(int ncid, int nvars, int varids[],
                   MPI_Offset *starts[], MPI_Offset *counts[],
                   void **buffers, MPI_Offset *bufcounts, 
                   MPI_Datatype *datatypes) {

  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems; 
  int *el_size, *nbytes;
  int total_nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype *ptype;
  MPI_Datatype buf_type;
  int isderived, *iscontig_of_ptypes;
  int i;
  MPI_Aint *displacement, a0, ai;
  int size;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;
  
  varp = (NC_var **)malloc(nvars*sizeof(NC_var *));  
  if (varp==NULL) printf("varp is NULL!!!\n"); 
  xbuf = (void **)malloc(nvars*sizeof(void *)); 
  if (xbuf==NULL) printf("xbuf is NULL!!!\n"); 
  cbuf = (void **)malloc(nvars*sizeof(void *)); 
  if (cbuf==NULL) printf("cbuf is NULL!!!\n"); 
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n"); 
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n"); 
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n"); 
  nbytes = malloc(nvars*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n"); 
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n"); 
  displacement = malloc(nvars*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n"); 
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n"); 
 
  for (i=0; i<nvars; i++){
  xbuf[i] = NULL;
  cbuf[i] = NULL; 
  varp[i] = ncmpii_NC_lookupvar(ncp, varids[i]);
  if(varp[i] == NULL)
    return NC_ENOTVAR;

  MPI_Type_size(datatypes[i], &size);
  status = ncmpii_dtype_decode(datatypes[i], &ptype[i], &el_size[i],
			       &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcounts[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (counts[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= counts[i][dim];
  }
  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;
  }
  /* set the mpi file view */
  
  status = set_vara_fileview_all(ncp, &(ncp->nciop->collective_fh), nvars, varp, starts, counts, 0);
  if(status != NC_NOERR)
    return status;

  total_nbytes = 0;
  for (i=0; i<nvars; i++) {/*loop for multi-variables*/
    if (!iscontig_of_ptypes[i]) {

      /* handling for derived datatype: pack into a contiguous buffer */

      cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );
      if (cbuf[i]==NULL) printf("cubf[%d] is NULL\n", i);
      status = ncmpii_data_repack((void *)buffers[i], bufcounts[i], datatypes[i], 
				cbuf[i], cnelems[i], ptype[i]);
      if (status != NC_NOERR)
        return status;
  
    } else {

      cbuf[i] = (void *)buffers[i];

    }

    /* assign or allocate MPI buffer */
    if ( need_convert(varp[i]->type, ptype[i]) ) {

      /* allocate new buffer */
      xbuf[i] = (void *)malloc(nbytes[i]);
      if (xbuf[i]==NULL) printf("xubf[%d] is NULL\n", i);

      /* automatic numeric datatype conversion */
   
      switch( varp[i]->type ) {
        case NC_BYTE:
           status = x_putn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_SHORT:
           status = x_putn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_INT:
           status = x_putn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_FLOAT:
           status = x_putn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_DOUBLE:
           status = x_putn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        default:
           break;
      } 

    } else  if ( need_swap(varp[i]->type, ptype[i]) ) {

      in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));
      xbuf[i] = (void *)cbuf[i];

    } else { 

      /* else, just assign contiguous buffer */
      xbuf[i] = (void *)cbuf[i];

    }
    total_nbytes = total_nbytes + nbytes[i];
    if (i==0){
	 displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else { 
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }
  }
 
  MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
  MPI_Type_commit(&buf_type);
  
  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf[0], 1, buf_type, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write_all");
        status = NC_EWRITE;
  }

  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  for (i=0; i<nvars; i++)
      if ( need_swap(varp[i]->type, ptype[i]) && buffers[i] == cbuf[i] && cbuf[i] == xbuf[i] )
          in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

    for (i=0; i<nvars; i++) {
        if (xbuf[i] != cbuf[i] && xbuf[i] != NULL) {
            free(xbuf[i]);
            xbuf[i] = NULL;
        }
        if (!iscontig_of_ptypes[i]) {
            if (cbuf[i] != buffers[i] && cbuf[i] != NULL) {
                free(cbuf[i]);
                cbuf[i] = NULL;
            }
        }
        if (status == NC_NOERR && IS_RECVAR(varp[i])) {
            /* update the number of records in header NC */
            MPI_Offset newnumrecs = starts[i][0] + counts[i][0];
            update_numrecs(ncp, newnumrecs);
        }
    }/*end for i*/

  free(varp);

  free(nbytes);
  free(el_size);

  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);

  free(xbuf);
  free(cbuf);

  MPI_Type_free(&buf_type);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_mvara_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts, 
                   MPI_Datatype *datatypes) {
  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems;
  int *el_size, *nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype *ptype;
  int isderived, *iscontig_of_ptypes;
  int i;
  int size;
  MPI_Aint *displacement, a0, ai;
  MPI_Datatype buf_type;
  

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  if(NC_indef(ncp))
    return NC_EINDEFINE;

  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
  varp = (NC_var **) malloc(nvars*sizeof(NC_var *));;
  if (varp==NULL) printf("varp is NULL!!!\n");
  xbuf = (void **)malloc(nvars*sizeof(void *));
  if (xbuf==NULL) printf("xbuf is NULL!!!\n");
  cbuf = (void **)malloc(nvars*sizeof(void *));
  if (cbuf==NULL) printf("cbuf is NULL!!!\n");
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n");
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n");
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n");
  nbytes = malloc(nvars*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n");
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n");
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n");
  displacement = malloc(nvars*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n");

  for (i=0; i<nvars; i++){
    xbuf[i] = NULL;
    cbuf[i] = NULL;

    varp[i] = ncmpii_NC_lookupvar(ncp, varids[i]);
    if(varp[i] == NULL)
     return NC_ENOTVAR;

    status = ncmpii_dtype_decode(datatypes[i], &ptype[i], &el_size[i],
                               &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
    if (status != NC_NOERR)
      return status;
    
    MPI_Type_size(datatypes[i], &size);
 
    if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcounts[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (counts[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= counts[i][dim];
  }

  if (nelems[i] != cnelems[i]) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems[i]>cnelems[i]) ? (nelems[i]=cnelems[i]) : (cnelems[i]=nelems[i]);
  }

  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;

  }
  /* set the mpi file view */

  status = set_vara_fileview_all(ncp, &(ncp->nciop->collective_fh), nvars, varp, starts, counts, 1);
//  status = set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp[i], start[i], count[i], 1);
  if(status != NC_NOERR)
    return status;

  for (i=0; i<nvars; i++){
  if (!iscontig_of_ptypes[i]) {

    /* account for derived datatype: allocate the contiguous buffer */

    cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );

  } else {

    cbuf[i] = (void *)buffers[i];

  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp[i]->type, ptype[i]) ) {

    /* allocate new buffer */
    xbuf[i] = (void *)malloc(nbytes[i]);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf[i] = (void *)cbuf[i];

  }
  if (i==0){
         displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else {
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }

  }

    MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
    MPI_Type_commit(&buf_type);

    mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf[0], 1, buf_type, &mpistatus);

//  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf[i], nbytes[i], MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        ncmpii_handle_error(rank, mpireturn, "MPI_File_read_all");
        status = NC_EREAD;
  }
  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  for (i=0; i<nvars; i++){
  if ( need_convert(varp[i]->type, ptype[i]) ) {

    /* automatic numeric datatype conversion */

    switch( varp[i]->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_INT:
         status = x_getn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp[i]->type, ptype[i]) ) {

    in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

  }

  if (!iscontig_of_ptypes[i]) {

    /* handling for derived datatype: unpack from the contiguous buffer */

    status = ncmpii_data_repack(cbuf[i], cnelems[i], ptype[i],
                                (void *)buffers[i], bufcounts[i], datatypes[i]);
    if (status != NC_NOERR)
      return status;

  }

  if (xbuf[i] != cbuf[i] && xbuf[i] != NULL)
    free(xbuf[i]);
  if (!iscontig_of_ptypes[i]) {
    if (cbuf[i] != buffers[i] && cbuf[i] != NULL)
      free(cbuf[i]);
    }
  } 
  free(varp);
  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_get_mvara_uchar_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_UNSIGNED_CHAR;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_schar_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_BYTE;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_text_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_CHAR;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_short_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_SHORT;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_int_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_INT;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_long_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_LONG;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_float_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_FLOAT;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_get_mvara_double_all(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
        datatypes[i] = MPI_DOUBLE;
  }
  ncmpi_get_mvara_all(ncid, nvars, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mget_vara_uchar_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_UNSIGNED_CHAR;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_schar_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_BYTE;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_text_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_CHAR;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_short_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_SHORT;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_int_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_INT;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_long_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_LONG;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_float_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_FLOAT;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

int
ncmpi_mget_vara_double_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts){
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int *varids;
  int i;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
 
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;
 
  for (i=0; i<ntimes; i++){
        datatypes[i] = MPI_DOUBLE;
        varids[i] = varid;
  }

  ncmpi_get_mvara_all(ncid, ntimes, varids,
                   starts, counts,
                   buffers, bufcounts,
                   datatypes);

  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }

  return (status);

}

static int
ncmpi_put_mvara_all_record(int ncid, int nvars, int varid[],
                   MPI_Offset *start[], MPI_Offset *count[],
                   void **buf, MPI_Offset *bufcount,
                   MPI_Datatype *datatype)
{

  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems;
  int *el_size, *nbytes;
  int total_nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype *ptype;
  MPI_Datatype buf_type;
  int isderived, *iscontig_of_ptypes;
  int i,j;
  MPI_Aint *displacement, a0, ai;
  int size;

  MPI_Barrier(MPI_COMM_WORLD);
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  comm = ncp->nciop->comm;

  if(NC_readonly(ncp))
    return NC_EPERM;

  if(NC_indef(ncp))
    return NC_EINDEFINE;

  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;

  varp = (NC_var **)malloc(nvars*sizeof(NC_var *));
  if (varp==NULL) printf("varp is NULL!!!\n");
  xbuf = (void **)malloc(nvars*sizeof(void *));
  if (xbuf==NULL) printf("xbuf is NULL!!!\n");
  cbuf = (void **)malloc(nvars*sizeof(void *));
  if (cbuf==NULL) printf("cbuf is NULL!!!\n");
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n");
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n");
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n");
  nbytes = malloc(nvars*count[0][0]*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n");
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n");
  displacement = malloc(nvars*count[0][0]*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n");
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n");

  for (i=0; i<nvars; i++){
  xbuf[i] = NULL;
  cbuf[i] = NULL;
  varp[i] = ncmpii_NC_lookupvar(ncp, varid[i]);
  if(varp[i] == NULL)
    return NC_ENOTVAR;
  MPI_Type_size(datatype[i], &size);
  status = ncmpii_dtype_decode(datatype[i], &ptype[i], &el_size[i],
                               &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcount[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (count[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= count[i][dim];
  }
 /* if (nelems[i] != cnelems[i]) {
 *     if (warning == NC_NOERR)
 *           warning = NC_EIOMISMATCH;
 *               (nelems[i]>cnelems[i]) ? (nelems[i]=cnelems[i]) : (cnelems[i]=nelems[i]);
 *                 }
 *                 */
  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;
  }
  /* set the mpi file view */

  status = set_vara_fileview_all(ncp, &(ncp->nciop->collective_fh), nvars, varp, start, count, 0);
  if(status != NC_NOERR)
    return status;

  total_nbytes = 0;
  for (i=0; i<nvars; i++) {/*loop for multi-variables*/
    if (!iscontig_of_ptypes[i]) {

      /* handling for derived datatype: pack into a contiguous buffer */

      cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );
      if (cbuf[i]==NULL) printf("cubf[%d] is NULL\n", i);
      status = ncmpii_data_repack((void *)buf[i], bufcount[i], datatype[i],
                                cbuf[i], cnelems[i], ptype[i]);
      if (status != NC_NOERR)
        return status;

    } else {

      cbuf[i] = (void *)buf[i];

    }

    /* assign or allocate MPI buffer */

    if ( need_convert(varp[i]->type, ptype[i]) ) {

      /* allocate new buffer */
      xbuf[i] = (void *)malloc(nbytes[i]);
      if (xbuf[i]==NULL) printf("xubf[%d] is NULL\n", i);

      /* automatic numeric datatype conversion */

      //printf("varp[i]->type:%d, ptype[i]:%d\n", varp[i]->type, ptype[i]);
      switch( varp[i]->type ) {
        case NC_BYTE:
           status = x_putn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_SHORT:
           status = x_putn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_INT:
           status = x_putn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_FLOAT:
           status = x_putn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_DOUBLE:
           status = x_putn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        default:
           break;
      }

    } else  if ( need_swap(varp[i]->type, ptype[i]) ) {

      in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));
      xbuf[i] = (void *)cbuf[i];

    } else {
      /* else, just assign contiguous buffer */
      xbuf[i] = (void *)cbuf[i];
    }
    total_nbytes = total_nbytes + nbytes[i];
    if (i==0){
         displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else {
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }
  }

  if (ncp->recsize <= varp[0]->len) {
    MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
  } else {
    for (i=0; i<nvars; i++){
        nbytes[i]=nbytes[i]/count[i][0];
    }
    for (j=1; j<count[0][0];j++){
      for (i=0; i<nvars; i++){
        nbytes[i+j*nvars]=nbytes[i];
        displacement[i+j*nvars] = displacement[i]+nbytes[i]*j;	
	}
    }
    MPI_Type_create_hindexed(nvars*count[0][0], nbytes, displacement, MPI_BYTE, &buf_type);
  }
  MPI_Type_commit(&buf_type);
  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf[0], 1, buf_type, &mpistatus);

  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        ncmpii_handle_error(rank, mpireturn, "MPI_File_write_all");
        status = NC_EWRITE;
  }
  /* reset fileview so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->collective_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  for (i=0; i<nvars; i++)
      if ( need_swap(varp[i]->type, ptype[i]) )
          in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

    for (i=0; i<nvars; i++){
        if (xbuf[i] != cbuf[i] && xbuf[i] != NULL) {
            free(xbuf[i]);
            xbuf[i]=NULL;
        }
        if (!iscontig_of_ptypes[i]) {
            if (cbuf[i] != buf[i] && cbuf[i] != NULL) {
                free(cbuf[i]);
                cbuf[i]=NULL;
            }
        }
        if (status == NC_NOERR && IS_RECVAR(varp[i])) {
            /* update the number of records in header NC */
            MPI_Offset newnumrecs = start[i][0] + count[i][0];
            update_numrecs(ncp, newnumrecs);
        }
    }/*end for i*/

  free(varp);

  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);

  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_mvara_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts,
                   MPI_Datatype *datatypes) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

 if (IS_RECVAR(varp))
  status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
 else 
  status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
 return (status);
}

static int
ncmpi_put_mvara_record(int ncid, int nvars, int varid[],
                   MPI_Offset *start[], MPI_Offset *count[],
                   void **buf, MPI_Offset *bufcount,
                   MPI_Datatype *datatype)
{

  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems;
  int *el_size, *nbytes;
  int total_nbytes;
  MPI_Status mpistatus;
//  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype *ptype;
  MPI_Datatype buf_type;
  int isderived, *iscontig_of_ptypes;
  int i,j;
  MPI_Aint *displacement, a0, ai;
  int size;

//  MPI_Barrier(MPI_COMM_WORLD);
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

//  comm = ncp->nciop->comm;

  if(NC_readonly(ncp))
    return NC_EPERM;

  if(NC_indef(ncp))
    return NC_EINDEFINE;

  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;

  varp = (NC_var **)malloc(nvars*sizeof(NC_var *));
  if (varp==NULL) printf("varp is NULL!!!\n");
  xbuf = (void **)malloc(nvars*sizeof(void *));
  if (xbuf==NULL) printf("xbuf is NULL!!!\n");
  cbuf = (void **)malloc(nvars*sizeof(void *));
  if (cbuf==NULL) printf("cbuf is NULL!!!\n");
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n");
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n");
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n");
  nbytes = malloc(nvars*count[0][0]*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n");
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n");
  displacement = malloc(nvars*count[0][0]*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n");
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n");

  for (i=0; i<nvars; i++){
  xbuf[i] = NULL;
  cbuf[i] = NULL;
  varp[i] = ncmpii_NC_lookupvar(ncp, varid[i]);
  if(varp[i] == NULL)
    return NC_ENOTVAR;
  MPI_Type_size(datatype[i], &size);
  status = ncmpii_dtype_decode(datatype[i], &ptype[i], &el_size[i],
                               &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcount[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (count[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= count[i][dim];
  }
 /* if (nelems[i] != cnelems[i]) {
 *     if (warning == NC_NOERR)
 *           warning = NC_EIOMISMATCH;
 *               (nelems[i]>cnelems[i]) ? (nelems[i]=cnelems[i]) : (cnelems[i]=nelems[i]);
 *                 }
 *                 */
  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;
  }
  /* set the mpi file view */

  status = set_vara_fileview_all(ncp, &(ncp->nciop->independent_fh), nvars, varp, start, count, 0);
//  status = set_vara_fileview(ncp, &(ncp->nciop->independent_fh), nvars, varp, start, count, 0);
  if(status != NC_NOERR)
    return status;

  total_nbytes = 0;
  for (i=0; i<nvars; i++) {/*loop for multi-variables*/
    if (!iscontig_of_ptypes[i]) {

      /* handling for derived datatype: pack into a contiguous buffer */

      cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );
      if (cbuf[i]==NULL) printf("cubf[%d] is NULL\n", i);
      status = ncmpii_data_repack((void *)buf[i], bufcount[i], datatype[i],
                                cbuf[i], cnelems[i], ptype[i]);
      if (status != NC_NOERR)
        return status;

    } else {

      cbuf[i] = (void *)buf[i];

    }

    /* assign or allocate MPI buffer */

    if ( need_convert(varp[i]->type, ptype[i]) ) {

      /* allocate new buffer */
      xbuf[i] = (void *)malloc(nbytes[i]);
      if (xbuf[i]==NULL) printf("xubf[%d] is NULL\n", i);

      /* automatic numeric datatype conversion */

      //printf("varp[i]->type:%d, ptype[i]:%d\n", varp[i]->type, ptype[i]);
      switch( varp[i]->type ) {
        case NC_BYTE:
           status = x_putn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_SHORT:
           status = x_putn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_INT:
           status = x_putn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_FLOAT:
           status = x_putn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_DOUBLE:
           status = x_putn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        default:
           break;
      }

    } else  if ( need_swap(varp[i]->type, ptype[i]) ) {

      in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));
      xbuf[i] = (void *)cbuf[i];

    } else {
      /* else, just assign contiguous buffer */
      xbuf[i] = (void *)cbuf[i];
    }
    total_nbytes = total_nbytes + nbytes[i];
    if (i==0){
         displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else {
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }
  }

  if (ncp->recsize <= varp[0]->len) {
    MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
  } else {
    for (i=0; i<nvars; i++){
        nbytes[i]=nbytes[i]/count[i][0];
    }
    for (j=1; j<count[0][0];j++){
      for (i=0; i<nvars; i++){
        nbytes[i+j*nvars]=nbytes[i];
        displacement[i+j*nvars] = displacement[i]+nbytes[i]*j;	
	}
    }
    MPI_Type_create_hindexed(nvars*count[0][0], nbytes, displacement, MPI_BYTE, &buf_type);
  }
  MPI_Type_commit(&buf_type);
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf[0], 1, buf_type, &mpistatus);

  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        ncmpii_handle_error(rank, mpireturn, "MPI_File_write");
        status = NC_EWRITE;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  for (i=0; i<nvars; i++)
      if ( need_swap(varp[i]->type, ptype[i]) )
          in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

  for (i=0; i<nvars; i++){
  if (xbuf[i] != cbuf[i] && xbuf[i] != NULL){
    free(xbuf[i]);
    xbuf[i]=NULL;
  }
  if (!iscontig_of_ptypes[i]) {
    if (cbuf[i] != buf[i] && cbuf[i] != NULL){
      free(cbuf[i]);
      cbuf[i]=NULL;
    }
  }


  if (status == NC_NOERR && IS_RECVAR(varp[i])) {
    MPI_Offset newnumrecs;
    newnumrecs = start[i][0] + count[i][0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
  }/*end for i*/

  free(varp);

  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);

  return ((warning != NC_NOERR) ? warning : status);
}

static int
ncmpi_put_mvara_nonrecord(int ncid, int nvars, int varids[],
                   MPI_Offset *starts[], MPI_Offset *counts[],
                   void **buffers, MPI_Offset *bufcounts, 
                   MPI_Datatype *datatypes) {

  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems; 
  int *el_size, *nbytes;
  int total_nbytes;
  MPI_Status mpistatus;
//  MPI_Comm comm;
  int mpireturn;
  MPI_Datatype *ptype;
  MPI_Datatype buf_type;
  int isderived, *iscontig_of_ptypes;
  int i;
  MPI_Aint *displacement, a0, ai;
  int size;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
//  comm = ncp->nciop->comm;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
  
  varp = (NC_var **)malloc(nvars*sizeof(NC_var *));  
  if (varp==NULL) printf("varp is NULL!!!\n"); 
  xbuf = (void **)malloc(nvars*sizeof(void *)); 
  if (xbuf==NULL) printf("xbuf is NULL!!!\n"); 
  cbuf = (void **)malloc(nvars*sizeof(void *)); 
  if (cbuf==NULL) printf("cbuf is NULL!!!\n"); 
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n"); 
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n"); 
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n"); 
  nbytes = malloc(nvars*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n"); 
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n"); 
  displacement = malloc(nvars*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n"); 
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n"); 
 
  for (i=0; i<nvars; i++){
  xbuf[i] = NULL;
  cbuf[i] = NULL; 
  varp[i] = ncmpii_NC_lookupvar(ncp, varids[i]);
  if(varp[i] == NULL)
    return NC_ENOTVAR;

  MPI_Type_size(datatypes[i], &size);
  status = ncmpii_dtype_decode(datatypes[i], &ptype[i], &el_size[i],
			       &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
  if (status != NC_NOERR)
    return status;

  if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcounts[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (counts[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= counts[i][dim];
  }
  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;
  }
  /* set the mpi file view */
  
  status = set_vara_fileview_all(ncp, &(ncp->nciop->independent_fh), nvars, varp, starts, counts, 0);
  if(status != NC_NOERR)
    return status;

  total_nbytes = 0;
  for (i=0; i<nvars; i++) {/*loop for multi-variables*/
    if (!iscontig_of_ptypes[i]) {

      /* handling for derived datatype: pack into a contiguous buffer */

      cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );
      if (cbuf[i]==NULL) printf("cubf[%d] is NULL\n", i);
      status = ncmpii_data_repack((void *)buffers[i], bufcounts[i], datatypes[i], 
				cbuf[i], cnelems[i], ptype[i]);
      if (status != NC_NOERR)
        return status;
  
    } else {

      cbuf[i] = (void *)buffers[i];

    }

    /* assign or allocate MPI buffer */
    if ( need_convert(varp[i]->type, ptype[i]) ) {

      /* allocate new buffer */
      xbuf[i] = (void *)malloc(nbytes[i]);
      if (xbuf[i]==NULL) printf("xubf[%d] is NULL\n", i);

      /* automatic numeric datatype conversion */
   
      switch( varp[i]->type ) {
        case NC_BYTE:
           status = x_putn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_SHORT:
           status = x_putn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_INT:
           status = x_putn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_FLOAT:
           status = x_putn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        case NC_DOUBLE:
           status = x_putn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
           break;
        default:
           break;
      } 

    } else  if ( need_swap(varp[i]->type, ptype[i]) ) {

      in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));
      xbuf[i] = (void *)cbuf[i];

    } else { 

      /* else, just assign contiguous buffer */
      xbuf[i] = (void *)cbuf[i];

    }
    total_nbytes = total_nbytes + nbytes[i];
    if (i==0){
	 displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else { 
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }
  }
 
  MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
  MPI_Type_commit(&buf_type);
  
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf[0], 1, buf_type, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
	ncmpii_handle_error(rank, mpireturn, "MPI_File_write");
        status = NC_EWRITE;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  for (i=0; i<nvars; i++)
      if ( need_swap(varp[i]->type, ptype[i]) && buffers[i] == cbuf[i] && cbuf[i] == xbuf[i] )
          in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

  for (i=0; i<nvars; i++){
  if (xbuf[i] != cbuf[i] && xbuf[i] != NULL){
    free(xbuf[i]);
    xbuf[i]=NULL;
  }
  if (!iscontig_of_ptypes[i]) {
    if (cbuf[i] != buffers[i] && cbuf[i] != NULL){
      free(cbuf[i]);
      cbuf[i]=NULL;
    }
  }

  if (status == NC_NOERR && IS_RECVAR(varp[i])) {
 
    /* update the number of records in NC
        and write it back to file header, if necessary
    */
    MPI_Offset newnumrecs;
    newnumrecs = starts[i][0] + counts[i][0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
  }/*end for i*/

  free(varp);

  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);

  return ((warning != NC_NOERR) ? warning : status);
}


int
ncmpi_put_mvara(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts,
                   MPI_Datatype *datatypes) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

 if (IS_RECVAR(varp))
  status = ncmpi_put_mvara_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
 else
  status = ncmpi_put_mvara_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
 return (status);
}


int
ncmpi_get_mvara(int ncid, int nvars, int *varids,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts, 
                   MPI_Datatype *datatypes) {
  NC_var **varp;
  NC *ncp;
  void **xbuf = NULL, **cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  MPI_Offset *nelems, *cnelems;
  int *el_size, *nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  MPI_Datatype *ptype;
  int isderived, *iscontig_of_ptypes;
  int i;
  int size;
  MPI_Aint *displacement, a0, ai;
  MPI_Datatype buf_type;
  

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  if(NC_indef(ncp))
    return NC_EINDEFINE;

  /* check to see that the desired mpi file handle is opened */

  status = check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
  varp = (NC_var **) malloc(nvars*sizeof(NC_var *));;
  if (varp==NULL) printf("varp is NULL!!!\n");
  xbuf = (void **)malloc(nvars*sizeof(void *));
  if (xbuf==NULL) printf("xbuf is NULL!!!\n");
  cbuf = (void **)malloc(nvars*sizeof(void *));
  if (cbuf==NULL) printf("cbuf is NULL!!!\n");
  nelems = malloc(nvars*sizeof(MPI_Offset));
  if (nelems==NULL) printf("nelems is NULL!!!\n");
  cnelems = malloc(nvars*sizeof(MPI_Offset));
  if (cnelems==NULL) printf("cnelems is NULL!!!\n");
  el_size = malloc(nvars*sizeof(int));
  if (el_size==NULL) printf("el_size is NULL!!!\n");
  nbytes = malloc(nvars*sizeof(int));
  if (nbytes==NULL) printf("nbytes is NULL!!!\n");
  ptype = malloc(nvars*sizeof(MPI_Datatype));
  if (ptype==NULL) printf("ptype is NULL!!!\n");
  iscontig_of_ptypes = malloc(nvars*sizeof(int));
  if (iscontig_of_ptypes==NULL) printf("iscontig_of_ptypes is NULL!!!\n");
  displacement = malloc(nvars*sizeof(MPI_Aint));
  if (displacement==NULL) printf("displacement is NULL!!!\n");

  for (i=0; i<nvars; i++){
    xbuf[i] = NULL;
    cbuf[i] = NULL;

    varp[i] = ncmpii_NC_lookupvar(ncp, varids[i]);
    if(varp[i] == NULL)
     return NC_ENOTVAR;

    status = ncmpii_dtype_decode(datatypes[i], &ptype[i], &el_size[i],
                               &cnelems[i], &isderived, &iscontig_of_ptypes[i]);
    if (status != NC_NOERR)
      return status;
    
    MPI_Type_size(datatypes[i], &size);
 
    if ( echar(varp[i]->type, ptype[i]) )
    return NC_ECHAR;

  cnelems[i] *= bufcounts[i];

  nelems[i] = 1;
  for (dim = 0; dim < varp[i]->ndims; dim++) {
    if (counts[i][dim] < 0)
      return NC_ENEGATIVECNT;
    nelems[i] *= counts[i][dim];
  }

  if (nelems[i] != cnelems[i]) {
    if (warning == NC_NOERR)
      warning = NC_EIOMISMATCH;
    (nelems[i]>cnelems[i]) ? (nelems[i]=cnelems[i]) : (cnelems[i]=nelems[i]);
  }

  nbytes[i] = varp[i]->xsz * nelems[i];
  if (nbytes[i] < 0)
    return NC_ENEGATIVECNT;

  }
  /* set the mpi file view */

  status = set_vara_fileview_all(ncp, &(ncp->nciop->independent_fh), nvars, varp, starts, counts, 1);
//  status = set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp[i], start[i], count[i], 1);
  if(status != NC_NOERR)
    return status;

  for (i=0; i<nvars; i++){
  if (!iscontig_of_ptypes[i]) {

    /* account for derived datatype: allocate the contiguous buffer */

    cbuf[i] = (void *)malloc( cnelems[i] * el_size[i] );

  } else {

    cbuf[i] = (void *)buffers[i];

  }

  /* assign or allocate MPI buffer */

  if ( need_convert(varp[i]->type, ptype[i]) ) {

    /* allocate new buffer */
    xbuf[i] = (void *)malloc(nbytes[i]);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf[i] = (void *)cbuf[i];

  }
  if (i==0){
         displacement[i] = 0;
         MPI_Get_address( xbuf[i], &a0 );
    } else {
         MPI_Get_address( xbuf[i], &ai );
         displacement[i] = (MPI_Aint) (ai-a0);
    }

  }

    MPI_Type_create_hindexed(nvars, nbytes, displacement, MPI_BYTE, &buf_type);
    MPI_Type_commit(&buf_type);

    mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf[0], 1, buf_type, &mpistatus);

//  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf[i], nbytes[i], MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        int rank;
        MPI_Comm_rank(ncp->nciop->comm, &rank);
        ncmpii_handle_error(rank, mpireturn, "MPI_File_read");
        status = NC_EREAD;
  }

  /* reset the file view so the entire file is visible again */
  MPI_File_set_view(ncp->nciop->independent_fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
 
  for (i=0; i<nvars; i++){
  if ( need_convert(varp[i]->type, ptype[i]) ) {

    /* automatic numeric datatype conversion */

    switch( varp[i]->type ) {
      case NC_BYTE:
         status = x_getn_schar(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_SHORT:
         status = x_getn_short(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_INT:
         status = x_getn_int(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_FLOAT:
         status = x_getn_float(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      case NC_DOUBLE:
         status = x_getn_double(xbuf[i], cbuf[i], cnelems[i], ptype[i]);
         break;
      default:
         break;
    }

  } else if ( need_swap(varp[i]->type, ptype[i]) ) {

    in_swapn(cbuf[i], nelems[i], ncmpix_len_nctype(varp[i]->type));

  }

  if (!iscontig_of_ptypes[i]) {

    /* handling for derived datatype: unpack from the contiguous buffer */

    status = ncmpii_data_repack(cbuf[i], cnelems[i], ptype[i],
                                (void *)buffers[i], bufcounts[i], datatypes[i]);
    if (status != NC_NOERR)
      return status;

  }

  if (xbuf[i] != cbuf[i] && xbuf[i] != NULL)
    free(xbuf[i]);
  if (!iscontig_of_ptypes[i]) {
    if (cbuf[i] != buffers[i] && cbuf[i] != NULL)
      free(cbuf[i]);
    }
  } 
  free(varp);
  free(nelems);
  free(cnelems);
  free(ptype);
  free(displacement);
  free(iscontig_of_ptypes);
 
  return ((warning != NC_NOERR) ? warning : status);
}

int
ncmpi_put_mvara_uchar_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_UNSIGNED_CHAR;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_schar_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] =  MPI_BYTE;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_text_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_CHAR;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_short_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_SHORT;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_int_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_INT;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_long_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_LONG;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_float_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_FLOAT;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_put_mvara_double_all(int ncid, int nvars, int varids[],
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varids[0]);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(nvars*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  for (i=0; i<nvars; i++){
	datatypes[i] = MPI_DOUBLE;	
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, nvars, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_uchar_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_UNSIGNED_CHAR;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_schar_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_BYTE;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_text_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_CHAR;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_short_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_SHORT;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_int_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_INT;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_long_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_LONG;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_float_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_FLOAT;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_mput_vara_double_all(int ncid, int ntimes, int varid,
                   MPI_Offset **starts, MPI_Offset **counts,
                   void **buffers, MPI_Offset *bufcounts) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
  MPI_Datatype *datatypes;
  int i;
  int *varids;
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  datatypes = (MPI_Datatype *)malloc(ntimes*sizeof(MPI_Datatype));
  if (datatypes == NULL) return NC_ENOMEM;
  
  varids = (int *)malloc(ntimes*sizeof(int));
  if (varids == NULL) return NC_ENOMEM;

  for (i=0; i<ntimes; i++){
	datatypes[i] = MPI_DOUBLE;	
	varids[i] = varid;
  }

  if (IS_RECVAR(varp))
    status = ncmpi_put_mvara_all_record(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  else
    status = ncmpi_put_mvara_all_nonrecord(ncid, ntimes, varids, starts, counts, buffers, bufcounts, datatypes);
  if (datatypes != NULL){
    free(datatypes);
    datatypes = NULL;
  }
  return (status);
}

int
ncmpi_iput_vara_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   MPI_Datatype datatype, NCMPI_Request *request) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;


  *request = ncmpii_new_NCMPI_Request(ncid, varid, varp->ndims,
		  start, count,
		  buf, bufcount, datatype);
  if (*request == NULL) printf("no memory buffer\n");

  return NC_NOERR;
}

int
ncmpi_iput_vara_uchar_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_UNSIGNED_CHAR, request);

  return status;
}

int
ncmpi_iput_vara_schar_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_BYTE, request);

  return status;
}

int
ncmpi_iput_vara_text_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_CHAR, request);

  return status;
}

int
ncmpi_iput_vara_short_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_SHORT, request);

  return status;
}

int
ncmpi_iput_vara_int_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_INT, request);

  return status;
}

int
ncmpi_iput_vara_long_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_LONG, request);

  return status;
}

int
ncmpi_iput_vara_float_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_FLOAT, request);

  return status;
}

int
ncmpi_iput_vara_double_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
   
  status = ncmpi_iput_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_DOUBLE, request);

  return status;
}

int
ncmpi_iget_vara_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   MPI_Datatype datatype, NCMPI_Request *request) {

  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;


  *request = ncmpii_new_NCMPI_Request(ncid, varid, varp->ndims,
		  start, count,
		  buf, bufcount, datatype);
  if (*request == NULL) printf("no memory buffer\n");
  return NC_NOERR;
}

int
ncmpi_iget_vara_uchar_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_UNSIGNED_CHAR, request);

  return status;
}

int
ncmpi_iget_vara_schar_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_BYTE, request);

  return status;
}

int
ncmpi_iget_vara_text_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_CHAR, request);

  return status;
}

int
ncmpi_iget_vara_short_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_SHORT, request);

  return status;
}

int
ncmpi_iget_vara_int_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_INT, request);

  return status;
}

int
ncmpi_iget_vara_float_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_FLOAT, request);

  return status;
}

int
ncmpi_iget_vara_double_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
                   const void *buf, MPI_Offset bufcount,
                   NCMPI_Request *request) {
  NC_var *varp;
  NC *ncp;
  int status = NC_NOERR;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = ncmpii_NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  status = ncmpi_iget_vara_all(ncid, varid, start, count, buf, bufcount,
                   MPI_DOUBLE, request);

  return status;
}

static int
ncmpi_coll_wait(NCMPI_Request request) {
  int ret;
  if (request->rw_flag == 1) {
    ret = ncmpi_put_vara_all(request->ncid, request->varid,
                   request->start, request->count,
                   request->buf, request->bufcount,
                   request->mpi_varatype);
  } else if ( request->rw_flag == 0) {
    ret = (ncmpi_get_vara_all(request->ncid, request->varid,
                   request->start, request->count,
                   request->buf, request->bufcount,
                   request->mpi_varatype));
  } else {
	  ret = NC_EFILE;
  }
  return ret;
}


int
ncmpi_coll_waitall(int count, NCMPI_Request array_of_requests[]) {
  int i,j;
  int ncid;
  int nvars;
  int *varids;
  MPI_Offset **starts;
  MPI_Offset **counts;
  void **buf;
  MPI_Offset *bufcount;
  MPI_Datatype *datatype;
  int ndim;

  ncid = array_of_requests[0]->ncid;
  nvars = count;
  varids = (int *)malloc(count*sizeof(int));
  starts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  counts = (MPI_Offset **)malloc(count*sizeof( MPI_Offset *));
  buf   = (void **)malloc(count*sizeof(void *));
  bufcount = (MPI_Offset *)malloc(count*sizeof(MPI_Offset));
  datatype = (MPI_Datatype *)malloc(count*sizeof(MPI_Datatype));
  
  for (i=0; i<count; i++){
        ndim = array_of_requests[i]->ndim;
	varids[i]=array_of_requests[i]->varid;
	starts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));	
	counts[i]=(MPI_Offset *)malloc(ndim*sizeof(MPI_Offset));	
	for (j=0; j<ndim; j++){
	   starts[i][j]=array_of_requests[i]->start[j];
	   counts[i][j]=array_of_requests[i]->count[j];
	}
        buf[i] = array_of_requests[i]->buf;
	bufcount[i] = array_of_requests[i]->bufcount;
	datatype[i] = array_of_requests[i]->mpi_varatype;
  }
  if (array_of_requests[0]->rw_flag == 1)
  ncmpi_put_mvara_all(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);

  if (array_of_requests[0]->rw_flag == 0)
  ncmpi_get_mvara_all(ncid, count, varids,
                   starts, counts,
                   buf, bufcount,
                   datatype);


  for (i=0; i<count; i++) {
	  free(starts[i]);
	  free(counts[i]);
	  ncmpii_free_NCMPI_Request(&(array_of_requests[i]));
  }
  free(starts);
  free(varids);
  free(counts);
  free(buf);
  free(bufcount);
  free(datatype);
  return NC_NOERR;
}

/* End non-blocking data access functions */
/* #################################################################### */
#endif
