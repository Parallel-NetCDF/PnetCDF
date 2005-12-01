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

#include "ncmpidtype.h"

/* Local prototypes */
static int length_of_mpitype(MPI_Datatype);

const char *
ncmpi_inq_libvers(void) {
  return "version = 1.0.0 of 27 July 2005";
}

/* Prototypes for functions used only in this file */
static int echar(nc_type nctype,MPI_Datatype mpitype);
static int need_convert(nc_type nctype,MPI_Datatype mpitype);
static int need_swap(nc_type nctype,MPI_Datatype mpitype);
static int x_putn_schar(void *xbuf, const void *buf, int nelems,
			MPI_Datatype datatype);
static int x_putn_short(void *xbuf, const void *buf, int nelems,
			MPI_Datatype datatype);
static int x_putn_int(void *xbuf, const void *buf, int nelems,
		      MPI_Datatype datatype);
static int x_putn_float(void *xbuf, const void *buf, int nelems,
			MPI_Datatype datatype);
static int x_putn_double(void *xbuf, const void *buf, int nelems,
			 MPI_Datatype datatype);
static int x_getn_schar(const void *xbuf, void *buf, int nelems,
			MPI_Datatype datatype);
static int x_getn_short(const void *xbuf, void *buf, int nelems,
			MPI_Datatype datatype);
static int x_getn_int(const void *xbuf, void *buf, int nelems,
		      MPI_Datatype datatype);
static int x_getn_float(const void *xbuf, void *buf, int nelems,
			MPI_Datatype datatype);
static int x_getn_double(const void *xbuf, void *buf, int nelems,
			 MPI_Datatype datatype);
static int set_var1_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp,
			     const MPI_Offset index[]);
static int set_var_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp);
static int set_vara_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp,
			     const MPI_Offset start[], const MPI_Offset count[],
			     int getnotput);
static int set_vars_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, 
			     const MPI_Offset start[], const MPI_Offset count[], 
			     const MPI_Offset stride[], int getnotput);
static int check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm,
		       int collective);

/* Begin Of Dataset Functions */

int 
ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp) {
  int status = NC_NOERR;
  size_t sizeof_off_t = 0;
  size_t chunksize=4098;	/* might be a good thing to hint later */
  NC *ncp;

  ncp = ncmpii_new_NC(&chunksize);
  if(ncp == NULL) 
    return NC_ENOMEM;

  assert(ncp->flags == 0);

  if (fIsSet(cmode, NC_64BIT_OFFSET)) {
	  /* unlike serial netcdf, we will not bother to support
	   * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
	   * serial netcdf has proven it's possible if datasets are small, but
	   * that's a hassle we don't want to worry about */
	  if (sizeof(off_t) != 8)
		  return NC_ESMALL;
	  fSet(cmode, NC_64BIT_OFFSET);
	  sizeof_off_t = 8;
  } else {
	  sizeof_off_t = 4;
  }
  assert(ncp->xsz = ncmpii_hdr_len_NC(ncp, sizeof_off_t));

  fSet(ncp->flags, NC_NOFILL);
  fSet(ncp->flags, NC_HSYNC);

  status = ncmpiio_create(comm, path, cmode, info, &ncp->nciop);  
  if(status != NC_NOERR) {
    ncmpii_free_NC(ncp);
    return status;
  }

  fSet(ncp->flags, NC_CREAT);

  if(fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
     /*
      * NC_SHARE implies sync up the number of records as well.
      * (File format version one.)
      * Note that other header changes are not shared
      * automatically.  Some sort of IPC (external to this package)
      * would be used to trigger a call to ncmpi_sync().
      */ 
    fSet(ncp->flags, NC_NSYNC);
  }

  ncmpii_add_to_NCList(ncp);
  *ncidp = ncp->nciop->fd;

  return status;
}

int
ncmpi_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, int *ncidp) {
  int status = NC_NOERR;
  NC *ncp;
  size_t chunksize=4098;	/* might be a good thing to hint later */
  
  ncp = ncmpii_new_NC(&chunksize);
  if(ncp == NULL)
    return NC_ENOMEM;

  status = ncmpiio_open(comm, path, omode, info, &ncp->nciop);
  if(status != NC_NOERR) {
    ncmpii_free_NC(ncp);
    return status;
  } 

  assert(ncp->flags == 0); 

  if(fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
    /*
     * NC_SHARE implies sync up the number of records as well.
     * (File format version one.)
     * Note that other header changes are not shared
     * automatically.  Some sort of IPC (external to this package)
     * would be used to trigger a call to ncmpi_sync().
     */ 
    fSet(ncp->flags, NC_NSYNC);
  }

  status = ncmpii_hdr_get_NC(ncp);
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
ncmpi_redef(int ncid) {
  int status;
  NC *ncp;
  int mynumrecs, numrecs;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR) 
    return status; 

  if(NC_readonly(ncp)) 
    return NC_EPERM;

  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* ensure exiting define mode always entering collective data mode */
  if(NC_indep(ncp))
    ncmpi_end_indep_data(ncid);

  if(fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
    /* read in from disk */
    status = ncmpii_read_NC(ncp);
    if(status != NC_NOERR)
      return status;
  } else {

    /* collect and set the max numrecs */

    mynumrecs = ncp->numrecs;
    MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_INT, MPI_MAX, ncp->nciop->comm);
    if (numrecs > ncp->numrecs) {
      ncp->numrecs = numrecs;
      set_NC_ndirty(ncp);
    }
  }

  ncp->old = ncmpii_dup_NC(ncp);
  if(ncp->old == NULL)
    return NC_ENOMEM;

  fSet(ncp->flags, NC_INDEF);

  return NC_NOERR;
}

int
ncmpi_begin_indep_data(int ncid) {
  int status = NC_NOERR;
  int mpireturn;
  NC *ncp;
  int rank;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if (status != NC_NOERR)
    return status;

  if (NC_indep(ncp))
    return NC_EINDEP;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

  if(!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
    mpireturn = MPI_File_sync(ncp->nciop->collective_fh);   /* collective */
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
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
  int status = NC_NOERR;
  int mpireturn;
  NC *ncp;
  int rank;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if (status != NC_NOERR)
    return status; 

  if (!NC_indep(ncp))
    return NC_ENOTINDEP;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  if(!NC_readonly(ncp) && NC_independentFhOpened(ncp->nciop)) {
    mpireturn = MPI_File_sync(ncp->nciop->independent_fh); /* independent */
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_sync error = %s\n", rank, errorString);
        MPI_Finalize();
        return NC_EFILE;
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
  if(status != NC_NOERR)
    return status;

  if(!NC_indef(ncp))
    return(NC_ENOTINDEFINE);

  return ncmpii_NC_enddef(ncp);
}

int
ncmpi_sync(int ncid) {
  int status = NC_NOERR;
  NC *ncp;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  if(NC_indef(ncp)) 
    return NC_EINDEFINE;

  if(NC_readonly(ncp))
    return ncmpii_read_NC(ncp);

  /* else, read/write */

  status = ncmpii_NC_sync(ncp);
  if(status != NC_NOERR)
    return status;

  return ncmpiio_sync(ncp->nciop);
}

int
ncmpi_abort(int ncid) {
 /*
  * In data mode, same as ncmpiio_close.
  * In define mode, descard new definition.
  * In create, remove the file.
  */
  int status;
  NC *ncp;
  int doUnlink = 0;

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
    status = ncmpii_NC_sync(ncp);
    if (status != NC_NOERR)
      return status;
  }

  (void) ncmpiio_close(ncp->nciop, doUnlink);
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
  if(status != NC_NOERR)
    return status;

  /* release NC object, close the file and write Dirty numrecs if necessary */

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
swapn(void *dst, const void *src, size_t nn, int xsize)
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

static int
x_putn_schar(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
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
x_putn_short(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
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
x_putn_int(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
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
x_putn_float(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
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
x_putn_double(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
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
x_getn_schar(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
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
x_getn_short(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
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
x_getn_int(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
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
x_getn_float(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
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
x_getn_double(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
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
      char errorString[512];
      int  errorStringLen;
      MPI_Error_string(mpireturn, errorString, &errorStringLen);
      printf("check_mpifh() calliing MPI_File_open error = %s\n", errorString);
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

/*
 * Check whether 'coord' values (indices) are valid for the variable.
 */
static int
NCcoordck(NC *ncp, const NC_var *varp, const MPI_Offset *coord)
{
        const MPI_Offset *ip;
        size_t *up;
 
        if(varp->ndims == 0)
                return NC_NOERR;        /* 'scalar' variable */
 
        if(IS_RECVAR(varp))
        {
                if(*coord > X_INT_MAX)
                        return NC_EINVALCOORDS; /* sanity check */
                if(NC_readonly(ncp) && *coord >= ncp->numrecs)
                {
                        if(!NC_doNsync(ncp))
                                return NC_EINVALCOORDS;
                        /* else */
                        {
                                /* Update from disk and check again */
                                const int status = ncmpii_read_numrecs(ncp);
                                if(status != NC_NOERR)
                                        return status;
                                if(*coord >= ncp->numrecs)
                                        return NC_EINVALCOORDS;
                        }
                }
                ip = coord + 1;
                up = varp->shape + 1;
        }
        else
        {
                ip = coord;
                up = varp->shape;
        }
 
        for(; ip < coord + varp->ndims; ip++, up++)
        {
                /* cast needed for braindead systems with signed size_t */
                if( *ip >= (unsigned long)*up )
                        return NC_EINVALCOORDS;
        }
 
        return NC_NOERR;                                                              }

/*
 * Check whether 'edges' are valid for the variable and 'start'
 */
/*ARGSUSED*/
static int
NCedgeck(const NC *ncp, const NC_var *varp,
         const MPI_Offset *start, const MPI_Offset *edges)
{
  const MPI_Offset *const end = start + varp->ndims;
  const size_t *shp = varp->shape;

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
    /* cast needed for braindead systems with signed size_t */
    if( *edges > (unsigned long)*shp || *start + *edges > (unsigned long)*shp)
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
  const size_t *shp = varp->shape; /* use size_t for now :( */

  if(varp->ndims == 0)
    return NC_NOERR;  /* 'scalar' variable */

  if(IS_RECVAR(varp))
  {
    if ( *stride == 0 || *stride >= X_INT_MAX)
      /* cast needed for braindead systems with signed size_t */
      return NC_ESTRIDE;

    start++;
    edges++;
    shp++;
    stride++;
  }

  for(; start < end; start++, edges++, shp++, stride++)
  {
    /* cast needed for braindead systems with signed size_t */
    if( (*edges > (unsigned long)*shp) || 
	(*edges > 0 && *start+1 + (*edges-1) * *stride > (unsigned long)*shp) ||
	(*edges == 0 && *start > (unsigned long)*shp) )
    {
      return(NC_EEDGE);
    }

    if ( *stride == 0 || *stride >= X_INT_MAX)
      /* cast needed for braindead systems with signed size_t */
      return NC_ESTRIDE;
  }

  return NC_NOERR;
}
#endif

static int
set_var1_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, const MPI_Offset index[]) {
  MPI_Offset offset;
  int status;
  int dim, ndims;
  int mpireturn;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  status = NCcoordck(ncp, varp, index);
  if (status != NC_NOERR)
    return status;

  offset = varp->begin;
 
  ndims = varp->ndims;

  if (ndims > 0) {

    if (IS_RECVAR(varp))
      offset += index[0] * ncp->recsize;
    else 
      offset += index[ndims-1] * varp->xsz;

    if (ndims > 1) {
      if (IS_RECVAR(varp))
        offset += index[ndims - 1] * varp->xsz;
      else
        offset += index[0] * varp->dsizes[1] * varp->xsz;

      for (dim = 1; dim < ndims - 1; dim++)
        offset += index[dim] * varp->dsizes[dim+1] * varp->xsz;
    }

  }

  mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, MPI_BYTE, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
        return NC_EFILE;
  }

  return NC_NOERR;
}

static int
set_var_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp) {
  MPI_Offset offset;
  int mpireturn;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  offset = varp->begin;

  if (!IS_RECVAR(varp)) { 
    /* Contiguous file view */
    mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, MPI_BYTE, "native", ncp->nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
        return NC_EFILE;
    }

  } else {
    /* Record variable, Strided file view */
    int  ndims;
    MPI_Datatype filetype;  
    MPI_Aint stride;
    int blocklen;

    ndims = varp->ndims;
    if (ndims > 1)
      blocklen = varp->dsizes[1] * varp->xsz;
    else
      blocklen = varp->xsz;

    stride = ncp->recsize;
    
    if (ncp->numrecs == 0)
	    return(NC_NOERR);

#if (MPI_VERSION < 2)
    MPI_Type_hvector(ncp->numrecs, blocklen, stride, MPI_BYTE, &filetype);
#else
    MPI_Type_create_hvector(ncp->numrecs, blocklen, stride, MPI_BYTE, &filetype);
#endif
    MPI_Type_commit(&filetype);

    mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, filetype, "native", ncp->nciop->mpiinfo);
    if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
        return NC_EFILE;
    }

    MPI_Type_free(&filetype); 
  }

  return NC_NOERR;
}

static int
set_vara_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, const MPI_Offset start[], const MPI_Offset count[], int getnotput) {

  MPI_Offset offset;
  int status;
  int dim, ndims;
  int *shape = NULL, *subcount = NULL, *substart = NULL; /* all in bytes */
  MPI_Datatype rectype;
  MPI_Datatype filetype;
  int mpireturn;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  offset = varp->begin;
  
  ndims = varp->ndims;

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

  if (ndims == 0) {

    /* scalar variable */
    filetype = MPI_BYTE;

  } else {

    /* if ndims == 0, all below pointers would be null */

    shape = (int *)malloc(sizeof(int) * ndims);
    subcount = (int *)malloc(sizeof(int) * ndims);
    substart = (int *)malloc(sizeof(int) * ndims);

    dim = 0;
    while (dim < ndims && count[dim] > 0) dim++;

    if (dim < ndims) {

      /* 0 size data */
      filetype = MPI_BYTE;

    } else {

      if (IS_RECVAR(varp)) {
        subcount[0] = count[0];
        substart[0] = 0;
        shape[0] = subcount[0];

        if (ncp->recsize <= varp->len) {
    
          /* the only record variable */

          if (varp->ndims == 1) {
            shape[0] *= varp->xsz;
	    subcount[0] *= varp->xsz;
          } else {
	    for (dim = 1; dim < ndims-1; dim++) {
              shape[dim] = varp->shape[dim];
              subcount[dim] = count[dim];
              substart[dim] = start[dim];
	    }
	    shape[dim] = varp->xsz * varp->shape[dim];
	    subcount[dim] = varp->xsz * count[dim];
	    substart[dim] = varp->xsz * start[dim];
          }
          offset += start[0] * ncp->recsize;

          MPI_Type_create_subarray(ndims, shape, subcount, substart, 
				    MPI_ORDER_C, MPI_BYTE, &filetype); 
    
          MPI_Type_commit(&filetype);
        } else {

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

          } else {
            for (dim = 1; dim < ndims-1; dim++) {
              shape[dim] = varp->shape[dim];
              subcount[dim] = count[dim];
              substart[dim] = start[dim];
            }
            shape[dim] = varp->xsz * varp->shape[dim];
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

      } else {

        /* non record variable */

        for (dim = 0; dim < ndims-1; dim++ ) {
          shape[dim] = varp->shape[dim];
          subcount[dim] = count[dim];
          substart[dim] = start[dim];
        }

        shape[dim] = varp->xsz * varp->shape[dim];
        subcount[dim] = varp->xsz * count[dim];
        substart[dim] = varp->xsz * start[dim];

        MPI_Type_create_subarray(ndims, shape, subcount, substart, 
		         MPI_ORDER_C, MPI_BYTE, &filetype); 

        MPI_Type_commit(&filetype);
      }

    }
  }

  mpireturn = MPI_File_set_view(*mpifh, offset, MPI_BYTE, 
		    filetype, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
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
  int rank;

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
        stride[dim] >= X_INT_MAX)
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
     (unsigned long)*start + (unsigned long)*count > NC_get_numrecs(ncp))
      return NC_EEDGE;
*/

  MPI_Comm_rank(ncp->nciop->comm, &rank);

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
                    *filetype, "native", ncp->nciop->mpiinfo);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_set_view error = %s\n", rank, errorString);
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
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

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

  /* set the mpi file view */
 
  status = set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if (status != NC_NOERR)
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf,  cbuf, 1, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes,
			     MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }
 
  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if ( status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

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

  /* set the mpi file view */
 
  status = set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if(status != NC_NOERR)
    return status; 

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

  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
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

    swapn(cbuf, xbuf, 1, ncmpix_len_nctype(varp->type));

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
ncmpi_get_var_all(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  /* set the mpi file view */
 
  status = set_var_fileview(ncp, &(ncp->nciop->collective_fh), varp);
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

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        status = NC_EREAD;
  }
 
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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
ncmpi_put_var(int ncid, int varid, const void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }
 
  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
ncmpi_get_var(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        status = NC_EREAD;
  }
 
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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
                   const void *buf, int bufcount, 
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
  MPI_Comm_rank(comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write_all error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
 
    /* update the number of records in NC
        and write it back to file header, if necessary
    */
    int newnumrecs, max_numrecs;
    newnumrecs = start[0] + count[0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
    if (NC_doNsync(ncp)) {
      MPI_Allreduce( &newnumrecs, &max_numrecs, 1, MPI_INT, MPI_MAX, comm );
 
      if (ncp->numrecs < max_numrecs) {
        ncp->numrecs = max_numrecs;
        if (rank == 0) {
          status = ncmpii_write_numrecs(ncp); /* call subroutine from nc.c */
          if(status != NC_NOERR)
            return status;
        }
      }
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}


int
ncmpi_get_vara_all(int ncid, int varid,
                   const MPI_Offset start[], const MPI_Offset count[],
		   void *buf, int bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  if ( need_convert(varp->type, ptype) ||  need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        status = NC_EREAD;
  }

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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
ncmpi_put_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL,  *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        status = NC_EREAD;
  }

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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
                   const void *buf, int bufcount, 
                   MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
  MPI_Comm_rank(comm, &rank);
 
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
 
  status = set_vars_fileview(ncp, &(ncp->nciop->collective_fh), 
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

  } else if ( need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write_all error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
 
    /* update the number of records in NC
        and write it back to file header, if necessary
    */
    int newnumrecs, max_numrecs;
    newnumrecs = start[0] + (count[0] - 1) * stride[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
    if (NC_doNsync(ncp)) {
      MPI_Allreduce( &newnumrecs, &max_numrecs, 1, MPI_INT, MPI_MAX, comm );
 
      if (ncp->numrecs < max_numrecs) {
        ncp->numrecs = max_numrecs;
        if (rank == 0) {
          status = ncmpii_write_numrecs(ncp); /* call subroutine from nc.c */
          if(status != NC_NOERR)
            return status;
        }
      }
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
}


int
ncmpi_get_vars_all(int ncid, int varid,
                   const MPI_Offset start[], 
		   const MPI_Offset count[],
                   const MPI_Offset stride[],
		   void *buf, int bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  if ( need_convert(varp->type, ptype) ||  need_swap(varp->type, ptype) ) {

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

  } else {

    /* else, just assign the contiguous buffer/user buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        status = NC_EREAD;
  }

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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }

  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        status = NC_EREAD;
  }
 
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

    swapn(cbuf, xbuf, nelems, ncmpix_len_nctype(varp->type));

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
		   const void *buf, int bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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

  status = ncmpi_put_vars_all(ncid, varid, start, count, stride,
                              cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    return status;

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
		   void *buf, int bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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

  status = ncmpi_get_vars_all(ncid, varid, start, count, stride, 
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
ncmpi_put_varm(int ncid, int varid,
	       const MPI_Offset start[],
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
	       const MPI_Offset imap[],
	       const void *buf, int bufcount,
	       MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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

  status = ncmpi_put_vars(ncid, varid, start, count, stride,
                          cbuf, cnelems, ptype);
  if (status != NC_NOERR)
    return status;

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
		   void *buf, int bufcount,
		   MPI_Datatype datatype) 
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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

/* #################################################################### */
/* Begin non-blocking data access functions */

void ncmpii_postwrite(void *xbuf, void *cbuf, void *buf) {
  if (xbuf != cbuf && xbuf != NULL)
    free(xbuf);
  if (cbuf != buf && cbuf != NULL)
    free(cbuf);
}

int ncmpii_postread(nc_type vartype,
		    void *xbuf, 
		    void *cbuf, 
		    int nelems, 
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

void ncmpii_postmwrite(void *cbuf, void *lbuf, int iscontig_of_ptypes) {
  if (!iscontig_of_ptypes && lbuf != NULL)
    free(lbuf);
  if (cbuf != NULL)
    free(cbuf);
}

int ncmpii_postmread(void *cbuf,
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

int 
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
ncmpi_wait(NCMPI_Request *request) {
  int mpireturn = MPI_SUCCESS;

  if (*request != NCMPI_REQUEST_NULL) {
    mpireturn = MPI_Wait(&((*request)->mpi_req), MPI_STATUS_IGNORE);
    ncmpii_postprocess(request);
  }

  if (mpireturn != MPI_SUCCESS) {
    return NC_EFILE;
  } else {
    return NC_NOERR;
  }
}

int 
ncmpi_waitall(int count, NCMPI_Request array_of_requests[]) {
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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
}

int
ncmpi_waitany(int count, NCMPI_Request array_of_requests[], int *index) {
  int i;
  int mpireturn = MPI_SUCCESS;
  MPI_Request *array_of_mpireqs;

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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

  array_of_mpireqs = (int *)malloc(count * sizeof(int));
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
               const void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

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

  /* set the mpi file view */
 
  status = set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if (status != NC_NOERR)
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf,  cbuf, 1, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, nbytes,
			      MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  if ( status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
               void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request) 
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

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

  /* set the mpi file view */
 
  status = set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if(status != NC_NOERR)
    return status; 

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

  mpireturn = MPI_File_iread(ncp->nciop->independent_fh, xbuf, 
			     nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
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
	       const void *buf, int bufcount, MPI_Datatype datatype,
	       NCMPI_Request *request) 
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }
 
  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;
 
  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
	       void *buf, int bufcount, MPI_Datatype datatype, 
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
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
ncmpi_iput_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               const void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL,  *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
    newnumrecs = start[0] + count[0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return ((warning != NC_NOERR) ? warning : status);
} 

int
ncmpi_iget_vara(int ncid, int varid,
               const MPI_Offset start[], const MPI_Offset count[],
               void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
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
ncmpi_iput_vars(int ncid, int varid,
               const MPI_Offset start[], 
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
               const void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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

    /* allocate new buffer */
    xbuf = (void *)malloc(nbytes);

    swapn(xbuf, cbuf, nelems, ncmpix_len_nctype(varp->type));

  } else {

    /* else, just assign contiguous buffer */
    xbuf = (void *)cbuf;

  }

  *request = (NCMPI_Request)malloc(sizeof(struct NCMPI_Req));

  mpireturn = MPI_File_iwrite(ncp->nciop->independent_fh, xbuf, 
			      nbytes, MPI_BYTE, &((*request)->mpi_req));
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        status = NC_EWRITE;
  }

  (*request)->reqtype = NCMPI_REQTYPE_WRITE;
  (*request)->xbuf = xbuf;
  (*request)->cbuf = cbuf;
  (*request)->buf = (void *)buf;
  (*request)->next_req = NCMPI_REQUEST_NULL;

  if (status == NC_NOERR && IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
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
               void *buf, int bufcount,
               MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  void *xbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int dim;
  int nelems, cnelems, el_size, nbytes;
  int mpireturn;
  int rank;
  MPI_Datatype ptype;
  int isderived, iscontig_of_ptypes;
 
  status = ncmpii_NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
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
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
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
ncmpi_iput_varm(int ncid, int varid,
	       const MPI_Offset start[],
	       const MPI_Offset count[],
	       const MPI_Offset stride[],
	       const MPI_Offset imap[],
	       const void *buf, int bufcount,
	       MPI_Datatype datatype,
	       NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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
		   void *buf, int bufcount,
		   MPI_Datatype datatype,
		   NCMPI_Request *request)
{
  NC_var *varp;
  NC *ncp;
  int ndims, dim;
  void *lbuf = NULL, *cbuf = NULL;
  int status = NC_NOERR, warning = NC_NOERR;
  int lnelems, cnelems, el_size;
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;
 
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
  int nelems;
 
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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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
  int nelems;

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

/* End non-blocking data access functions */
/* #################################################################### */
