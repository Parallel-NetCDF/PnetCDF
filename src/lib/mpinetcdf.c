/*********************************************************************************
 * 
 * This file is written by Northwestern University and Argonne National Laboratory
 *
 ********************************************************************************/

#include <mpi.h>
#include "nc.h"
#include "ncx.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>

const char *
ncmpi_inq_libvers(void) {
  return "version = 3.5.0 of Aug 30 2002 13:00:00 $";
}

/* to be updated */ 
const char *
ncmpi_strerror(int err) {
  return nc_strerror(err);
}

/* Begin Of Dataset Functions */

int 
ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp) {
  int status = NC_NOERR;
  NC *ncp;

  ncp = new_NC(NULL);
  if(ncp == NULL) 
    return NC_ENOMEM;

  assert(ncp->xsz = hdr_len_NC(ncp));
  assert(ncp->flags == 0);

  fSet(ncp->flags, NC_NOFILL);

  status = ncio_create(comm, path, cmode, info, &ncp->nciop);  
  if(status != NC_NOERR) {
    free_NC(ncp);
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

  add_to_NCList(ncp);
  *ncidp = ncp->nciop->fd;

  return status;
}

int
ncmpi_open(MPI_Comm comm, const char *path, int omode, MPI_Info info, int *ncidp) {
  int status = NC_NOERR;
  NC *ncp;
  
  ncp = new_NC(NULL);
  if(ncp == NULL)
    return NC_ENOMEM;

  status = ncio_open(comm, path, omode, info, &ncp->nciop);
  if(status != NC_NOERR) {
    free_NC(ncp);
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

  status = hdr_get_NC(ncp);
  if (status != NC_NOERR) {
    free_NC(ncp);
    return status;
  }

  add_to_NCList(ncp);

  *ncidp = ncp->nciop->fd;
 
  return status;
}

int
ncmpi_redef(int ncid) {
  int status;
  NC *ncp;
  int mynumrecs, numrecs;

  status = NC_check_id(ncid, &ncp);
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
    status = read_NC(ncp);
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

  ncp->old = dup_NC(ncp);
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

  status = NC_check_id(ncid, &ncp);
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
 
  status = NC_check_id(ncid, &ncp);
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

  status = NC_check_id(ncid, &ncp); 
  if(status != NC_NOERR)
    return status;

  if(!NC_indef(ncp))
    return(NC_ENOTINDEFINE);

  return NC_enddef(ncp);
}

int
ncmpi_sync(int ncid) {
  int status = NC_NOERR;
  NC *ncp;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  if(NC_indef(ncp)) 
    return NC_EINDEFINE;

  if(NC_readonly(ncp))
    return read_NC(ncp);

  /* else, read/write */

  status = NC_sync(ncp);
  if(status != NC_NOERR)
    return status;

  return ncio_sync(ncp->nciop);
}

int
ncmpi_abort(int ncid) {
 /*
  * In data mode, same as ncio_close.
  * In define mode, descard new definition.
  * In create, remove the file.
  */
  int status;
  NC *ncp;
  int doUnlink = 0;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  doUnlink = NC_IsNew(ncp);

  if (ncp->old != NULL) {
    /* a plain redef, not a create */
    assert(!NC_IsNew(ncp));
    assert(fIsSet(ncp->flags, NC_INDEF));
    free_NC(ncp->old);
    ncp->old = NULL;
    fClr(ncp->flags, NC_INDEF);
  } 
  else if (!NC_readonly(ncp) && !NC_indef(ncp)) {
    /* data mode, write */
    status = NC_sync(ncp);
    if (status != NC_NOERR)
      return status;
  }

  (void) ncio_close(ncp->nciop, doUnlink);
  ncp->nciop = NULL;

  del_from_NCList(ncp);

  free_NC(ncp);

  return NC_NOERR;
}

int
ncmpi_close(int ncid) {
  int status = NC_NOERR;
  NC *ncp;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  /* release NC object, close the file and write Dirty numrecs if necessary */

  return NC_close(ncp);
}

/* End Of Dataset Functions */

/* Begin Of Define Mode Functions */

int
ncmpi_def_dim(int ncid, const char *name, size_t len, int *idp) {
  return nc_def_dim(ncid, name, len, idp);
}

int
ncmpi_def_var(int ncid, const char *name, nc_type xtype, 
              int ndims, const int *dimidsp, int *varidp) {
  return nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp);
}

int 
ncmpi_rename_dim(int ncid, int dimid, const char *name) {
  return nc_rename_dim(ncid, dimid, name);
}

int 
ncmpi_rename_var(int ncid, int varid, const char *name) {
  return nc_rename_var(ncid, varid, name);
}

/* End Of Define Mode Functions */

/* Begin Of Inquiry Functions */

int 
ncmpi_inq(int ncid, int *ndimsp, int *nvarsp,
          int *ngattsp, int *unlimdimidp) {
  return nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);
}

int 
ncmpi_inq_ndims(int ncid, int *ndimsp) {
  return nc_inq_ndims(ncid, ndimsp);
}

int
ncmpi_inq_nvars(int ncid, int *nvarsp) {
  return nc_inq_nvars(ncid, nvarsp); 
}

int 
ncmpi_inq_natts(int ncid, int *ngattsp) {
  return nc_inq_natts(ncid, ngattsp);
} 

int
ncmpi_inq_unlimdim(int ncid, int *unlimdimidp) {
  return nc_inq_unlimdim(ncid, unlimdimidp);
} 

int 
ncmpi_inq_dimid(int ncid, const char *name, int *idp) {
  return nc_inq_dimid(ncid, name, idp);
} 

int 
ncmpi_inq_dim(int ncid, int dimid, char *name, size_t *lenp) {
  return nc_inq_dim(ncid, dimid, name, lenp);
} 

int 
ncmpi_inq_dimname(int ncid, int dimid, char *name) {
  return nc_inq_dimname(ncid, dimid, name);
} 

int 
ncmpi_inq_dimlen(int ncid, int dimid, size_t *lenp) {
  return nc_inq_dimlen(ncid, dimid, lenp);
} 

int 
ncmpi_inq_var(int ncid, int varid, char *name, 
              nc_type *xtypep, int *ndimsp, int *dimidsp,
              int *nattsp) {
  return nc_inq_var(ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp);
} 

int 
ncmpi_inq_varid(int ncid, const char *name, int *varidp) {
  return nc_inq_varid(ncid, name, varidp);
} 

int 
ncmpi_inq_varname(int ncid, int varid, char *name) {
  return nc_inq_varname(ncid, varid, name);
} 

int 
ncmpi_inq_vartype(int ncid, int varid, nc_type *xtypep) {
  return nc_inq_vartype(ncid, varid, xtypep);
} 

int 
ncmpi_inq_varndims(int ncid, int varid, int *ndimsp) {
  return nc_inq_varndims(ncid, varid, ndimsp);
} 

int 
ncmpi_inq_vardimid(int ncid, int varid, int *dimidsp) {
  return nc_inq_vardimid(ncid, varid, dimidsp);
} 

int
ncmpi_inq_varnatts(int ncid, int varid, int *nattsp) {
  return nc_inq_varnatts(ncid, varid, nattsp);
} 

/* End Of Inquiry Functions */

/* Begin Of Attribute Functions */

int 
ncmpi_inq_att(int ncid, int varid, const char *name,
              nc_type *xtypep, size_t *lenp) {
  return nc_inq_att(ncid, varid, name, xtypep, lenp);
}

int 
ncmpi_inq_attid(int ncid, int varid, const char *name, int *idp) {
  return nc_inq_attid(ncid, varid, name, idp);
}

int 
ncmpi_inq_atttype(int ncid, int varid, const char *name,
                  nc_type *xtypep) {
  return nc_inq_atttype(ncid, varid, name, xtypep);
}

int 
ncmpi_inq_attlen(int ncid, int varid, const char *name,
                 size_t *lenp) {
  return nc_inq_attlen(ncid, varid, name, lenp);
}

int 
ncmpi_inq_attname(int ncid, int varid, int attnum, char *name) {
  return nc_inq_attname(ncid, varid, attnum, name);
}

int 
ncmpi_copy_att(int ncid_in, int varid_in, const char *name,
               int ncid_out, int varid_out) {
  return nc_copy_att(ncid_in, varid_in, name, ncid_out, varid_out);
}

int 
ncmpi_rename_att(int ncid, int varid, const char *name,
                 const char *newname) {
  return nc_rename_att(ncid, varid, name, newname);
}

int 
ncmpi_del_att(int ncid, int varid, const char *name) {
  return nc_del_att(ncid, varid, name);
}

int 
ncmpi_put_att_text(int ncid, int varid, const char *name, size_t len,
                   const char *op) {
  return nc_put_att_text(ncid, varid, name, len, op);
}

int 
ncmpi_get_att_text(int ncid, int varid, const char *name, char *ip) {
  return nc_get_att_text(ncid, varid, name, ip);
}

int 
ncmpi_put_att_uchar(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const unsigned char *op) {
  return nc_put_att_uchar(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_uchar(int ncid, int varid, const char *name,
                    unsigned char *ip) {
  return nc_get_att_uchar(ncid, varid, name, ip);
}

int 
ncmpi_put_att_schar(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const signed char *op) {
  return nc_put_att_schar(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_schar(int ncid, int varid, const char *name,
                    signed char *ip) {
  return nc_get_att_schar(ncid, varid, name, ip);
}

int 
ncmpi_put_att_short(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const short *op) {
  return nc_put_att_short(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_short(int ncid, int varid, const char *name, short *ip) {
  return nc_get_att_short(ncid, varid, name, ip);
}

int 
ncmpi_put_att_int(int ncid, int varid, const char *name,
                  nc_type xtype, size_t len, const int *op) {
  return nc_put_att_int(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_int(int ncid, int varid, const char *name, int *ip) {
  return nc_get_att_int(ncid, varid, name, ip);
}

int 
ncmpi_put_att_long(int ncid, int varid, const char *name,
                   nc_type xtype, size_t len, const long *op) {
  return nc_put_att_long(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_long(int ncid, int varid, const char *name, long *ip) {
  return nc_get_att_long(ncid, varid, name, ip);
}

int 
ncmpi_put_att_float(int ncid, int varid, const char *name,
                    nc_type xtype, size_t len, const float *op) {
  return nc_put_att_float(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_float(int ncid, int varid, const char *name, float *ip) {
  return nc_get_att_float(ncid, varid, name, ip);
}

int 
ncmpi_put_att_double(int ncid, int varid, const char *name,
                     nc_type xtype, size_t len, const double *op) {
  return nc_put_att_double(ncid, varid, name, xtype, len, op);
}

int 
ncmpi_get_att_double(int ncid, int varid, const char *name, 
                     double *ip) {
  return nc_get_att_double(ncid, varid, name, ip);
}

/* End Of Attribute Functions */

/* Begin {put,get}_att */

/* to be updated */
#if 0
int
ncmpi_put_att_text(int ncid, int varid, const char *name,
        size_t len, const char *op) {
  int status = NC_NOERR;
/*
  int rank;
  MPI_Comm_rank(comm, &rank);
*/
  status = nc_put_att_text(ncid, varid, name, len, op);

  return status;
}
#endif

/* End {put,get}_att */

/* Begin {put,get}_var */

/*
 *  MAPPING:  MPI DATATYPE   <--->   NETCDF DATATYPE
 *		MPI_BYTE		NC_BYTE
 *		MPI_CHAR		NC_CHAR
 *		MPI_SHORT		NC_SHORT
 *		MPI_INT			NC_INT
 *		MPI_FLOAT		NC_FLOAT
 *		MPI_DOUBLE		NC_DOUBLE
 *
 *
 *  Assume: MPI_Datatype and nc_type are both enumerable types
 */

int
length_of_mpitype(MPI_Datatype datatype) {
  switch(datatype) {
    case MPI_BYTE:
    case MPI_CHAR:
	return ((int)sizeof(char));
    case MPI_SHORT:
	return (int)(sizeof(short));
    case MPI_INT:
	return((int)sizeof(int)); 
    case MPI_FLOAT:
	return((int)sizeof(float));
    case MPI_DOUBLE:
	return((int)sizeof(double));
  }

  return -1;
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

void
x_putn_short(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch (datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_put_short_short(xp, (const short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_put_short_int(xp, (const int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_put_short_float(xp, (const float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_put_short_double(xp, (const double *)data);
        return;
  }
} 

void
x_putn_int(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch (datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_put_int_short(xp, (const short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_put_int_int(xp, (const int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_put_int_float(xp, (const float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_put_int_double(xp, (const double *)data);
        return;
  }
} 

void
x_putn_float(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch (datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_put_float_short(xp, (const short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_put_float_int(xp, (const int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_put_float_float(xp, (const float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_put_float_double(xp, (const double *)data);
        return;
  }
} 

void
x_putn_double(void *xbuf, const void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf; 
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);  

  switch (datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_put_double_short(xp, (const short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_put_double_int(xp, (const int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_put_double_float(xp, (const float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc) 
          ncx_put_double_double(xp, (const double *)data);
        return;
  }
} 

void
x_getn_short(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;

  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch(datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_get_short_short(xp, (short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_get_short_int(xp, (int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_get_short_float(xp, (float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, data += datainc)
          ncx_get_short_double(xp, (double *)data);
        return;
  }
} 

void 
x_getn_int(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch(datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_get_int_short(xp, (short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_get_int_int(xp, (int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_get_int_float(xp, (float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, data += datainc)
          ncx_get_int_double(xp, (double *)data);
        return;
  }
} 

void
x_getn_float(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
 
  switch(datatype) {
    case MPI_SHORT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_get_float_short(xp, (short *)data);
        return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_get_float_int(xp, (int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_get_float_float(xp, (float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, data += datainc)
          ncx_get_float_double(xp, (double *)data);
        return;
  }
}

void
x_getn_double(const void *xbuf, void *buf, int nelems, MPI_Datatype datatype) {
  char *xp, *data;
  int datainc;
 
  xp = (char *) xbuf;
  data = (char *) buf;
  datainc = length_of_mpitype(datatype);
                                                                                    
  switch(datatype) {
    case MPI_SHORT:
	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc) 
	  ncx_get_double_short(xp, (short *)data);
	return;
    case MPI_INT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_get_double_int(xp, (int *)data);
        return;
    case MPI_FLOAT:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_get_double_float(xp, (float *)data);
        return;
    case MPI_DOUBLE:
        for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, data += datainc)
          ncx_get_double_double(xp, (double *)data);
        return;
  }
}

int
NC_check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm, int collective) {

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
      printf("NC_check_mpifh() calliing MPI_File_open error = %s\n", errorString);
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
int
NCcoordck(NC *ncp, const NC_var *varp, const size_t *coord)
{
        const size_t *ip;
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
                                const int status = read_numrecs(ncp);
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
                if((unsigned long) *ip >= (unsigned long) *up )
                        return NC_EINVALCOORDS;
        }
 
        return NC_NOERR;                                                              }

int
NC_set_var1_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, const size_t index[]) {
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

int
NC_set_var_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp) {
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

int
NC_set_vara_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, const size_t start[], const size_t count[]) {

  MPI_Offset offset;
  int status;
  int dim, ndims;
  int *shape, *subcount, *substart; /* all in bytes */
  size_t *end;
  MPI_Datatype rectype;
  MPI_Datatype filetype;
  int mpireturn;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  offset = varp->begin;
  
  ndims = varp->ndims;

  if (ndims > 0) {

    /* if ndims == 0, all below pointers would be null */

    end = (size_t *)malloc(sizeof(size_t) * ndims);
    shape = (int *)malloc(sizeof(int) * ndims);
    subcount = (int *)malloc(sizeof(int) * ndims);
    substart = (int *)malloc(sizeof(int) * ndims);

    for (dim = 0; dim < ndims; dim++)
      end[dim] = start[dim] + count[dim] - 1;
  }

  status = NCcoordck(ncp, varp, end);
  if (status != NC_NOERR)
    return status;
  
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
	MPI_Type_vector(subcount[0], varp->xsz, ncp->recsize,
			MPI_BYTE, &filetype);
#else
	MPI_Type_create_vector(subcount[0], varp->xsz, ncp->recsize,
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

    if (ndims > 0) {
      shape[dim] = varp->xsz * varp->shape[dim];
      subcount[dim] = varp->xsz * count[dim];
      substart[dim] = varp->xsz * start[dim];

      MPI_Type_create_subarray(ndims, shape, subcount, substart, 
	 		     MPI_ORDER_C, MPI_BYTE, &filetype); 

      MPI_Type_commit(&filetype);
    } else {
      /* scalar variable */
      filetype = MPI_BYTE;
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
    MPI_Type_free(&filetype);

    free(end);
    free(shape);
    free(subcount);
    free(substart);
  }

  return NC_NOERR;
}

int
NC_set_vars_fileview(NC* ncp, MPI_File *mpifh, NC_var* varp, 
		     const size_t start[], const size_t count[], 
                     const size_t stride[]) {
  MPI_Offset offset;
  int status;
  int mpireturn;
  int dim, ndims;
  MPI_Datatype *subtypes, *filetype;
  size_t *blocklens, *blockstride, *blockcount, *end;
  int rank;

  MPI_Comm_rank(ncp->nciop->comm, &rank);

  offset = varp->begin;
  
  ndims = varp->ndims;

  if (ndims == 0) {

    /* scalar variable */

    filetype = subtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype));
    *filetype = MPI_BYTE;
  
  } else {  
  
    subtypes = (MPI_Datatype *)malloc((ndims+1) * sizeof(MPI_Datatype));
    filetype = subtypes;
    subtypes[ndims] = MPI_BYTE;
  
    end = (size_t *) malloc(ndims * sizeof(size_t));
    for (dim = 0; dim < ndims; dim++)
      end[dim] = start[dim] + (count[dim] - 1) * stride[dim];
    status = NCcoordck(ncp, varp, end);
    if (status != NC_NOERR)
      return status; 
  
    blocklens = (size_t *) malloc(ndims * sizeof(size_t));
    blockstride = (size_t *) malloc(ndims * sizeof(size_t));
    blockcount = (size_t *) malloc(ndims * sizeof(size_t));
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
    for (dim=0; dim < ndims; dim++) 
      MPI_Type_free(subtypes + dim);

    free(blocklens);
    free(blockstride);
    free(blockcount);
  }

  free(subtypes);

  return NC_NOERR;
}

int
ncmpi_put_var1(int ncid, int varid,
               const size_t index[],
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int words_bigendian = 0;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if(status != NC_NOERR)
    return status; 

  nbytes = varp->xsz;

  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;	
  }
  else
  { 
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
   
    /* automatic datatype conversion */
   
    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, 1, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, 1, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, 1, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, 1, datatype);
           break;
        default:
           break;
      }
   
    } else if (!words_bigendian) {
   
      swapn(xbuf,  buf, 1, ncx_len_nctype(varp->type));
   
    }
  }
 
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

 
  if (xbuf != buf)
    free(xbuf);
 
  if (IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
    newnumrecs = index[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return status;
}

int
ncmpi_get_var1(int ncid, int varid,
               const size_t index[],
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;
 
#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);

  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  /* set the mpi file view */
 
  status = NC_set_var1_fileview(ncp, &(ncp->nciop->independent_fh), varp, index);
  if(status != NC_NOERR)
    return status; 

  nbytes = varp->xsz;

  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR ) 
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else 
  {  
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
 
  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        return NC_EREAD;
  }

 
  /* automatic datatype conversion */
 
  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
 
      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, 1, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, 1, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, 1, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, 1, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
 
    } else if (!words_bigendian) {
   
      swapn(buf, xbuf, 1, ncx_len_nctype(varp->type));
      free(xbuf);
 
    }
  }
 
  return status;
}

int
ncmpi_get_var_all(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = bufcount/length_of_mpitype(datatype);
  nbytes = nelems * varp->xsz;

  /* set the mpi file view */
 
  status = NC_set_var_fileview(ncp, &(ncp->nciop->collective_fh), varp);
  if(status != NC_NOERR)
    return status;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
 
  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        return NC_EREAD;
  }

 
  /* automatic datatype conversion */
 
  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
 
      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
  
    } else if (!words_bigendian) {
   
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);
   
    }
  }
 
  return status;
}

int
ncmpi_put_var(int ncid, int varid, const void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int nelems, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int words_bigendian = 0;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 
 
  nelems = bufcount/length_of_mpitype(datatype);
  nbytes = nelems * varp->xsz;

  /* set the mpi file view */
 
  status = NC_set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
  if(status != NC_NOERR)
    return status;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
   
    /* automatic datatype conversion */
   
    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
 
    } else if (!words_bigendian) {
 
      swapn(xbuf,  buf, nelems, ncx_len_nctype(varp->type));
 
    }
  }
 
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

 
  if (xbuf != buf)
    free(xbuf);
 
  if (IS_RECVAR(varp)) {
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
 
  return status;
}

int
ncmpi_get_var(int ncid, int varid, void *buf, int bufcount, MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = bufcount/length_of_mpitype(datatype);
  nbytes = nelems * varp->xsz;
 
  /* set the mpi file view */
 
  status = NC_set_var_fileview(ncp, &(ncp->nciop->independent_fh), varp);
  if(status != NC_NOERR)
    return status;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
 
  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        return NC_EREAD;
  }

 
  /* automatic datatype conversion */
 
  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
 
      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
 
    } else if (!words_bigendian) {
 
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);
 
    }
  }

  return status; 
} 

int
ncmpi_put_vara_all(int ncid, int varid,
                   const size_t start[], const size_t count[],
                   const void *buf, int bufcount, 
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  int rank;
  int words_bigendian = 0;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
  MPI_Comm_rank(comm, &rank);
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp, start, count);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;

  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);

    /* automatic datatype conversion */

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }

    } else if (!words_bigendian) { 
  
      swapn(xbuf,  buf, nelems, ncx_len_nctype(varp->type));
  
    }
  }

  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write_all error = %s\n", rank, errorString);
        return NC_EWRITE;
  }


  if (xbuf != buf)
    free(xbuf);
 
  if (IS_RECVAR(varp)) {
 
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
          status = write_numrecs(ncp); /* call subroutine from nc.c */
          if(status != NC_NOERR)
            return status;
        }
      }
    }
  }
 
  return status;
}


int
ncmpi_get_vara_all(int ncid, int varid,
                   const size_t start[], const size_t count[],
		   void *buf, int bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vara_fileview(ncp, &(ncp->nciop->collective_fh), varp, start, count);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems; 

  /* assign or allocate MPI read buffer */ 

  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian && 
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }

  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        return NC_EREAD;
  }


  /* automatic datatype conversion */

  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {

      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);

    } else if (!words_bigendian) {
  
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);

    }
  }
 
  return status;
}

int
ncmpi_put_vara(int ncid, int varid,
               const size_t start[], const size_t count[],
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int words_bigendian = 0;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp, start, count);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
 
    /* automatic datatype conversion */ 
   
    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
 
    } else if (!words_bigendian) {
 
      swapn(xbuf,  buf, nelems, ncx_len_nctype(varp->type));
 
    }
  }
 
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

  if (xbuf != buf)
    free(xbuf);

  if (IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
    newnumrecs = start[0] + count[0];
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }
 
  return status;
} 

int
ncmpi_get_vara(int ncid, int varid,
               const size_t start[], const size_t count[],
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vara_fileview(ncp, &(ncp->nciop->independent_fh), varp, start, count);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
 
  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        return NC_EREAD;
  }
 
  /* automatic datatype conversion */
 
  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
 
      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
 
    } else if (!words_bigendian) {
 
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);
 
    }
  }
 
  return status;
} 

int
ncmpi_put_vars_all(int ncid, int varid,
                   const size_t start[], 
		   const size_t count[],
		   const size_t stride[],
                   const void *buf, int bufcount, 
                   MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  MPI_Status mpistatus;
  MPI_Comm comm;
  int mpireturn;
  int rank;
  int words_bigendian = 0;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  comm = ncp->nciop->comm;
  MPI_Comm_rank(comm, &rank);
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->collective_fh), comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vars_fileview(ncp, &(ncp->nciop->collective_fh), 
                                varp, start, count, stride);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;

  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
 
    /* automatic datatype conversion */
    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
 
    } else if (!words_bigendian) {
 
      swapn(xbuf,  buf, nelems, ncx_len_nctype(varp->type));
 
    }
  }
 
  mpireturn = MPI_File_write_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write_all error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

 
  if (xbuf != buf)
    free(xbuf);
 
  if (IS_RECVAR(varp)) {
 
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
          status = write_numrecs(ncp); /* call subroutine from nc.c */
          if(status != NC_NOERR)
            return status;
        }
      }
    }
  }
 
  return status;
}


int
ncmpi_get_vars_all(int ncid, int varid,
                   const size_t start[], 
		   const size_t count[],
                   const size_t stride[],
		   void *buf, int bufcount,
                   MPI_Datatype datatype) {

  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->collective_fh), ncp->nciop->comm, 1);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vars_fileview(ncp, &(ncp->nciop->collective_fh), 
				varp, start, count, stride);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems; 

  /* assign or allocate MPI read buffer */ 

  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian && 
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
  
  mpireturn = MPI_File_read_all(ncp->nciop->collective_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read_all error = %s\n", rank, errorString);
        return NC_EREAD;
  }


  /* automatic datatype conversion */

  if ( ncx_len_nctype(varp->type) != 1 ) { 

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {

      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
  
    } else if (!words_bigendian) {
  
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);
  
    }
  }
 
  return status;
}

int
ncmpi_put_vars(int ncid, int varid,
               const size_t start[], 
	       const size_t count[],
	       const size_t stride[],
               const void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  MPI_Status mpistatus;
  int mpireturn;
  int words_bigendian = 0;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_readonly(ncp))
    return NC_EPERM;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
			        varp, start, count, stride);
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = (void *)buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
 
    /* automatic datatype conversion */ 
 
    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
      switch( varp->type ) {
        case NC_SHORT:
           x_putn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_putn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_putn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_putn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
   
    } else if (!words_bigendian) {
 
      swapn(xbuf,  buf, nelems, ncx_len_nctype(varp->type));
 
    }
  }
 
  mpireturn = MPI_File_write(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_write error = %s\n", rank, errorString);
        return NC_EWRITE;
  }

 
  if (xbuf != buf)
    free(xbuf);

  if (IS_RECVAR(varp)) {
    /* update the number of records in NC */
 
    int newnumrecs;
    newnumrecs = start[0] + (count[0] - 1) * stride[0] + 1;
    if (ncp->numrecs < newnumrecs) {
      ncp->numrecs = newnumrecs;
      set_NC_ndirty(ncp);
    }
  }

  return status;
} 

int
ncmpi_get_vars(int ncid, int varid,
               const size_t start[], 
	       const size_t count[],
               const size_t stride[],
               void *buf, int bufcount,
               MPI_Datatype datatype) {
  NC_var *varp;
  NC *ncp;
  void *xbuf;
  int status;
  int dim;
  int nelems, nbytes;
  int words_bigendian = 0;
  MPI_Status mpistatus;
  int mpireturn;
  int rank;

#if WORDS_BIGENDIAN
  words_bigendian = 1;
#endif
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  if(NC_indef(ncp))
    return NC_EINDEFINE;
 
  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  /* check to see that the desired mpi file handle is opened */
 
  status = NC_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  /* set the mpi file view */
 
  status = NC_set_vars_fileview(ncp, &(ncp->nciop->independent_fh),
				varp, start, count, stride); 
  if(status != NC_NOERR)
    return status;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = varp->xsz * nelems;
 
  /* assign or allocate MPI read buffer */
 
  if ( varp->type == NC_BYTE || varp->type == NC_CHAR )
  {
    assert( length_of_mpitype(datatype) == 1 );
    xbuf = (void *)buf;
  }
  else
  {
    if ( words_bigendian &&
        ncx_len_nctype(varp->type) == length_of_mpitype(datatype) )
      /* Just assign MPI read buffer */
      xbuf = buf;
    else
      /* else, allocate new buffer */
      xbuf = (void *)malloc(nbytes);
  }
 
  mpireturn = MPI_File_read(ncp->nciop->independent_fh, xbuf, nbytes, MPI_BYTE, &mpistatus);
  if (mpireturn != MPI_SUCCESS) {
        char errorString[512];
        int  errorStringLen;
        MPI_Error_string(mpireturn, errorString, &errorStringLen);
        printf("%2d: MPI_File_read error = %s\n", rank, errorString);
        return NC_EREAD;
  }

 
  /* automatic datatype conversion */
 
  if ( ncx_len_nctype(varp->type) != 1 ) {

    if ( ncx_len_nctype(varp->type) != length_of_mpitype(datatype) ) {
 
      switch( varp->type ) {
        case NC_SHORT:
           x_getn_short(xbuf, buf, nelems, datatype);
           break;
        case NC_INT:
           x_getn_int(xbuf, buf, nelems, datatype);
           break;
        case NC_FLOAT:
           x_getn_float(xbuf, buf, nelems, datatype);
           break;
        case NC_DOUBLE:
           x_getn_double(xbuf, buf, nelems, datatype);
           break;
        default:
           break;
      }
      free(xbuf);
   
    } else if (!words_bigendian) {
 
      swapn(buf, xbuf, nelems, ncx_len_nctype(varp->type));
      free(xbuf);
 
    }
  }
 
  return status;
} 

int
ncmpi_put_var1_text(int ncid, int varid,
                     const size_t index[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, sizeof(char), MPI_CHAR);
}


int
ncmpi_put_var1_short(int ncid, int varid,
                     const size_t index[],
		     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_put_var1(ncid, varid, index, 
                        (const void *)op, sizeof(short), MPI_SHORT); 
}

int
ncmpi_put_var1_int(int ncid, int varid,
                   const size_t index[],
                   const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, sizeof(int), MPI_INT);
}

int
ncmpi_put_var1_float(int ncid, int varid,
                     const size_t index[],
                     const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, sizeof(float), MPI_FLOAT); 
}
 
int
ncmpi_put_var1_double(int ncid, int varid,
                      const size_t index[],
                      const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_put_var1(ncid, varid, index,
                        (const void *)op, sizeof(double), MPI_DOUBLE);
}

int
ncmpi_get_var1_text(int ncid, int varid,
                     const size_t index[],
                     char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, sizeof(char), MPI_CHAR);
}

int
ncmpi_get_var1_short(int ncid, int varid,
                     const size_t index[],
                     short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, sizeof(short), MPI_SHORT); 
}
 
int
ncmpi_get_var1_int(int ncid, int varid,
                   const size_t index[],
                   int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, sizeof(int), MPI_INT);  
}
 
int
ncmpi_get_var1_float(int ncid, int varid,
                     const size_t index[],
                     float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, sizeof(float), MPI_FLOAT);  
}
 
int
ncmpi_get_var1_double(int ncid, int varid,
                      const size_t index[],
                      double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  return ncmpi_get_var1(ncid, varid, index,
                        (void *)ip, sizeof(double), MPI_DOUBLE);  
} 

int
ncmpi_put_var_text(int ncid, int varid, const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];

  nbytes = nelems * sizeof(char);

  return ncmpi_put_var(ncid, varid, (const void *)op, nbytes, MPI_CHAR);
}

int
ncmpi_put_var_short(int ncid, int varid, const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  ndims = varp->ndims;

  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp)) 
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];

  nbytes = nelems * sizeof(short);
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nbytes, MPI_SHORT);
}

int
ncmpi_put_var_int(int ncid, int varid, const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(int);
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nbytes, MPI_INT);
} 

int
ncmpi_put_var_float(int ncid, int varid, const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(float);
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nbytes, MPI_FLOAT);
} 

int
ncmpi_put_var_double(int ncid, int varid, const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(double);
 
  return ncmpi_put_var(ncid, varid, (const void *)op, nbytes, MPI_DOUBLE);
} 

int
ncmpi_get_var_text(int ncid, int varid, char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];

  nbytes = nelems * sizeof(char);

  return ncmpi_get_var(ncid, varid, (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_var_short(int ncid, int varid, short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(short);
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nbytes, MPI_SHORT);
}

int
ncmpi_get_var_int(int ncid, int varid, int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(int);
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nbytes, MPI_INT);
} 

int
ncmpi_get_var_float(int ncid, int varid, float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(float);
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nbytes, MPI_FLOAT);
} 

int
ncmpi_get_var_double(int ncid, int varid, double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(double);
 
  return ncmpi_get_var(ncid, varid, (void *)ip, nbytes, MPI_DOUBLE);
} 

int
ncmpi_get_var_text_all(int ncid, int varid, char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  ndims = varp->ndims;

  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];

  nbytes = nelems * sizeof(char);

  return ncmpi_get_var_all(ncid, varid, (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_var_short_all(int ncid, int varid, short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(short);
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nbytes, MPI_SHORT);
}

int
ncmpi_get_var_int_all(int ncid, int varid, int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(int);
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nbytes, MPI_INT);
} 

int
ncmpi_get_var_float_all(int ncid, int varid, float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(float);
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nbytes, MPI_FLOAT);
} 

int
ncmpi_get_var_double_all(int ncid, int varid, double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int ndims;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  ndims = varp->ndims;
 
  if (ndims > 1)
    nelems = varp->dsizes[1];
  else
    nelems = 1;
  if (IS_RECVAR(varp))
    nelems *= ncp->numrecs;
  else
    nelems *= varp->shape[0];
 
  nbytes = nelems * sizeof(double);
 
  return ncmpi_get_var_all(ncid, varid, (void *)ip, nbytes, MPI_DOUBLE);
} 

int
ncmpi_put_vara_text_all(int ncid, int varid,
                         const size_t start[], const size_t count[],
                         const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nbytes, MPI_CHAR);
}

int
ncmpi_put_vara_text(int ncid, int varid,
                     const size_t start[], const size_t count[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nbytes, MPI_CHAR);
}

int
ncmpi_put_vara_short_all(int ncid, int varid,
                         const size_t start[], const size_t count[],
                         const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;
 
  return ncmpi_put_vara_all(ncid, varid, start, count,
                            (const void *)op, nbytes, MPI_SHORT);
}

int
ncmpi_put_vara_short(int ncid, int varid,
                     const size_t start[], const size_t count[],
                     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nbytes, MPI_SHORT);
} 

int
ncmpi_put_vara_int_all(int ncid, int varid, 
                       const size_t start[], const size_t count[], 
		       const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;                                                                  

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR; 

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;

  return ncmpi_put_vara_all(ncid, varid, start, count, 
			    (const void *)op, nbytes, MPI_INT);
}

int
ncmpi_put_vara_int(int ncid, int varid, 
		const size_t start[], const size_t count[], 
		const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nbytes, MPI_INT);
}

int
ncmpi_put_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;

  return ncmpi_put_vara_all(ncid, varid, start, count, 
			    (const void *)op, nbytes, MPI_FLOAT);
}

int
ncmpi_put_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nbytes, MPI_FLOAT);
}

int
ncmpi_put_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_put_vara_all(ncid, varid, start, count, 
                            (const void *)op, nbytes, MPI_DOUBLE);
}

int
ncmpi_put_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                const double *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;
 
  return ncmpi_put_vara(ncid, varid, start, count,
                        (const void *)op, nbytes, MPI_DOUBLE);
}

int
ncmpi_get_vara_text_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_vara_text(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    char *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_vara_short_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    short *ip) {
 
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;
 
  return ncmpi_get_vara_all(ncid, varid, start, count,
                            (void *)ip, nbytes, MPI_SHORT);
}

int
ncmpi_get_vara_short(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    short *ip) {
 
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nbytes, MPI_SHORT);
} 

int
ncmpi_get_vara_int_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    int *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;
  
  return ncmpi_get_vara_all(ncid, varid, start, count, 
                            (void *)ip, nbytes, MPI_INT);
}

int
ncmpi_get_vara_int(int ncid, int varid,
                const size_t start[], const size_t count[],
                int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nbytes, MPI_INT);
}

int
ncmpi_get_vara_float_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    float *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;
  
  return ncmpi_get_vara_all(ncid, varid, start, count, 
			    (void *)ip, nbytes, MPI_FLOAT); 
}

int
ncmpi_get_vara_float(int ncid, int varid,
                const size_t start[], const size_t count[],
                float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nbytes, MPI_FLOAT);
}

int
ncmpi_get_vara_double_all(int ncid, int varid,
                    const size_t start[], const size_t count[],
                    double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_get_vara_all(ncid, varid, start, count, 
			    (void *)ip, nbytes, MPI_DOUBLE); 
}

int
ncmpi_get_vara_double(int ncid, int varid,
                const size_t start[], const size_t count[],
                double *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;
 
  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
 
  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;
 
  return ncmpi_get_vara(ncid, varid, start, count,
                        (void *)ip, nbytes, MPI_DOUBLE);
}

int
ncmpi_put_vars_text_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nbytes, MPI_CHAR);
}

int
ncmpi_put_vars_text(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     const char *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nbytes, MPI_CHAR);
}

int
ncmpi_put_vars_short_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nbytes, MPI_SHORT);
}

int
ncmpi_put_vars_short(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     const short *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nbytes, MPI_SHORT);
}

int
ncmpi_put_vars_int_all(int ncid, int varid,
                       const size_t start[],
                       const size_t count[],
                       const size_t stride[],
                       const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nbytes, MPI_INT);
}

int
ncmpi_put_vars_int(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   const int *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;
  
  nelems = 1; 
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;
    
  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nbytes, MPI_INT);
} 

int
ncmpi_put_vars_float_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nbytes, MPI_FLOAT);
}

int
ncmpi_put_vars_float(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     const float *op) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nbytes, MPI_FLOAT);
}

int
ncmpi_put_vars_double_all(int ncid, int varid,
                          const size_t start[], 
                          const size_t count[],
                          const size_t stride[],
                          const double *op) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_put_vars_all(ncid, varid, start, count, stride,
                            (const void *)op, nbytes, MPI_DOUBLE);

}

int
ncmpi_put_vars_double(int ncid, int varid,
                      const size_t start[],
                      const size_t count[],
                      const size_t stride[],
                      const double *op) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_put_vars(ncid, varid, start, count, stride,
                        (const void *)op, nbytes, MPI_DOUBLE);

}

int
ncmpi_get_vars_text_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_vars_text(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     char *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(char) * nelems;

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nbytes, MPI_CHAR);
}

int
ncmpi_get_vars_short_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nbytes, MPI_SHORT);
}

int
ncmpi_get_vars_short(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     short *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(short) * nelems;

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nbytes, MPI_SHORT);
}

int
ncmpi_get_vars_int_all(int ncid, int varid,
                       const size_t start[],
                       const size_t count[],
                       const size_t stride[],
                       int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nbytes, MPI_INT);
}

int
ncmpi_get_vars_int(int ncid, int varid,
                   const size_t start[],
                   const size_t count[],
                   const size_t stride[],
                   int *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(int) * nelems;

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nbytes, MPI_INT);
}

int
ncmpi_get_vars_float_all(int ncid, int varid,
                         const size_t start[],
                         const size_t count[],
                         const size_t stride[],
                         float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
                            (void *)ip, nbytes, MPI_FLOAT);
}

int
ncmpi_get_vars_float(int ncid, int varid,
                     const size_t start[],
                     const size_t count[],
                     const size_t stride[],
                     float *ip) {
  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(float) * nelems;

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nbytes, MPI_FLOAT);
}

int
ncmpi_get_vars_double_all(int ncid, int varid,
                          const size_t start[], 
			  const size_t count[],
			  const size_t stride[],
                          double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;
 
  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_get_vars_all(ncid, varid, start, count, stride,
			    (void *)ip, nbytes, MPI_DOUBLE); 
}

int
ncmpi_get_vars_double(int ncid, int varid,
                      const size_t start[],
                      const size_t count[],
                      const size_t stride[],
                      double *ip) {

  NC_var *varp;
  NC *ncp;
  int status;
  int dim;
  int nelems, nbytes;

  status = NC_check_id(ncid, &ncp);
  if(status != NC_NOERR)
    return status;

  varp = NC_lookupvar(ncp, varid);
  if(varp == NULL)
    return NC_ENOTVAR;

  nelems = 1;
  for (dim = 0; dim < varp->ndims; dim++)
    nelems *= count[dim];
  nbytes = (int)sizeof(double) * nelems;

  return ncmpi_get_vars(ncid, varid, start, count, stride,
                        (void *)ip, nbytes, MPI_DOUBLE);
}

/* End {put,get}_var */
