/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include "nc.h"
#include "rnd.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include <strings.h> /* bzero() */
#include <assert.h>
#include "ncx.h"
#include "macro.h"

/* list of open netcdf's */
static NC *NClist = NULL;

/* This is the default create format for ncmpi_create and nc__create. */
static int default_create_format = NC_FORMAT_CLASSIC;

/* These have to do with version numbers. */
#define MAGIC_NUM_LEN 4
#define VER_CLASSIC 1
#define VER_64BIT_OFFSET 2
#define VER_HDF5 3
#define VER_64BIT_DATA 5




/* Prototypes for functions used only in this file */
static int move_data_r(NC *ncp, NC *old);
static int move_recs_r(NC *ncp, NC *old);
static int move_vars_r(NC *ncp, NC *old);
static int write_NC(NC *ncp);
static int NC_begins(NC *ncp, MPI_Offset h_minfree, MPI_Offset v_align,
                     MPI_Offset v_minfree, MPI_Offset r_align);
static int NC_check_header(MPI_Comm comm, void *buf, MPI_Offset nn, NC *ncp);
static int NC_check_def(MPI_Comm comm, void *buf, MPI_Offset nn);

#if 0
static int enddef(NC *ncp);
static int nc_sync(int ncid);
static int nc_set_fill(int ncid, int fillmode, int *old_mode_ptr);
#endif

void
ncmpii_add_to_NCList(NC *ncp)
{
        assert(ncp != NULL);

        ncp->prev = NULL;
        if(NClist != NULL)
                NClist->prev = ncp;
        ncp->next = NClist;
        NClist = ncp;
}

void
ncmpii_del_from_NCList(NC *ncp)
{
        assert(ncp != NULL);

        if(NClist == ncp)
        {
                assert(ncp->prev == NULL);
                NClist = ncp->next;
        }
        else
        {
                assert(ncp->prev != NULL);
                ncp->prev->next = ncp->next;
        }

        if(ncp->next != NULL)
                ncp->next->prev = ncp->prev;

        ncp->next = NULL;
        ncp->prev = NULL;
}

/*
 * Check the data set definitions across all processes by
 * comparing the header buffer streams of all processes.
 */
static int
NC_check_header(MPI_Comm comm, void *buf, MPI_Offset hsz, NC *ncp) {
    int rank, errcheck, status=NC_NOERR, errflag, compare;
    MPI_Offset hsz_0;
    void *cmpbuf;
    bufferinfo gbp;

    MPI_Comm_rank(comm, &rank);

    /* process 0 broadcasts header size */
    if (rank == 0) hsz_0 = hsz;
    MPI_Bcast(&hsz_0, 1, MPI_LONG_LONG_INT, 0, comm);

    if (rank == 0)
        cmpbuf = buf;
    else
        cmpbuf = (void*) NCI_Malloc(hsz_0);

    /* process 0 broadcasts its header */
    MPI_Bcast(cmpbuf, hsz_0, MPI_BYTE, 0, comm);

    compare = 0;
    if (hsz != hsz_0)  /* hsz may be different from hsz_0 */
        compare = 1;
    else if (rank > 0)
        compare = memcmp(buf, cmpbuf, hsz);

    /* Use LOR because memcmp() may return negative value */
    MPI_Allreduce(&compare, &errcheck, 1, MPI_INT, MPI_LOR, comm);
    if (errcheck == 0) {
        if (rank > 0)
            NCI_Free(cmpbuf);
        return NC_NOERR;
    }
    /* now part of the header is not consistent across all processes */

    gbp.nciop  = ncp->nciop;
    gbp.offset = 0;    /* read from start of the file */
    gbp.size   = hsz_0;
    gbp.index  = 0;
    gbp.pos    = gbp.base = cmpbuf;

    status = ncmpii_hdr_check_NC(&gbp, ncp);

    errflag = (status == NC_NOERR) ? 0 : 1;
    MPI_Allreduce(&errflag, &errcheck, 1, MPI_INT, MPI_MAX, comm);
    if (errcheck > 0) {
        status = (status != NC_NOERR) ? status : NC_EMULTIDEFINE;
    }
    else if (hsz != hsz_0) {
        /* TODO !!!! need to replace the local NC object with root's.
           The local header size is different from the root's and the
           header consistency check returns only a warning. We must
           overwrite local header, NC object pointed by ncp, with
           the root's header, current in cmpbuf */
    }
    if (rank > 0) NCI_Free(cmpbuf);

    return NC_NOERR;
}


/* 'defined but not used': seems like a useful function though. why did we
 * write it?  should we be using it? */

static int
NC_check_def(MPI_Comm comm, void *buf, MPI_Offset nn) {
  int rank;
  int errcheck;
  MPI_Offset compare = 0;
  void *cmpbuf;
  MPI_Offset max_size;

  MPI_Comm_rank(comm, &rank);

  if (rank == 0)
    max_size = nn;
  MPI_Bcast(&max_size, 1, MPI_LONG_LONG_INT, 0, comm);

  compare = max_size - nn;

  MPI_Allreduce(&compare, &errcheck, 1, MPI_LONG_LONG_INT, MPI_LOR, comm);

  if (errcheck)
    return NC_EMULTIDEFINE;

  if (rank == 0)
    cmpbuf = buf;
  else
    cmpbuf = (void *)NCI_Malloc(nn);

  MPI_Bcast(cmpbuf, nn, MPI_BYTE, 0, comm);

  if (rank != 0) {
    compare = memcmp(buf, cmpbuf, nn);
    NCI_Free(cmpbuf);
  }

  MPI_Allreduce(&compare, &errcheck, 1, MPI_LONG_LONG_INT, MPI_LOR, comm);

  if (errcheck){
    return NC_EMULTIDEFINE;
  }else{
    return NC_NOERR;
  }
}

int
ncmpii_NC_check_id(int ncid, NC **ncpp)
{
        NC *ncp;

        if(ncid >= 0)
        {
                for(ncp = NClist; ncp != NULL; ncp = ncp->next)
                {
                        if(ncp->nciop->fd == ncid)
                        {
                                *ncpp = ncp;
                                return NC_NOERR; /* normal return */
                        }
                }
        }

        /* else, not found */
        return NC_EBADID;
}


/* static */
void
ncmpii_free_NC(NC *ncp)
{
        if(ncp == NULL)
                return;
        ncmpii_free_NC_dimarrayV(&ncp->dims);
        ncmpii_free_NC_attrarrayV(&ncp->attrs);
        ncmpii_free_NC_vararrayV(&ncp->vars);
        NCI_Free(ncp);
}


/* static */
NC *
ncmpii_new_NC(const MPI_Offset *chunkp)
{
        NC *ncp;

        ncp = (NC *) NCI_Malloc(sizeof(NC));
        if(ncp == NULL)
                return NULL;
        bzero(ncp, sizeof(NC));

        if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
                ncp->xsz = MIN_NC_XSZ+28; /*Kgao*/
        } else {
                ncp->xsz = MIN_NC_XSZ; 
        }

        assert(ncp->xsz == ncmpii_hdr_len_NC(ncp, 0)); 
        
        ncp->chunk = (chunkp != NULL) ? *chunkp : NC_SIZEHINT_DEFAULT;

        return ncp;
}

/* This function sets a default create flag that will be logically
   or'd to whatever flags are passed into nc_create for all future
   calls to nc_create.
   Valid default create flags are NC_64BIT_OFFSET, NC_CLOBBER,
   NC_LOCK, NC_SHARE. */
int
ncmpi_set_default_format(int format, int *old_formatp)
{
    /* Return existing format if desired. */
    if (old_formatp)
      *old_formatp = default_create_format;

    /* Make sure only valid format is set. */
    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT && format != NC_FORMAT_64BIT_DATA  ){
      return NC_EINVAL;
    }
    default_create_format = format;
    return NC_NOERR;
}

#if 0
/* returns a value suituable for a create flag.  Will return one or more of the
 * following values ORed together:
 * NC_64BIT_OFFSET, NC_CLOBBER, NC_LOCK, NC_SHARE */
static int
ncmpii_get_default_format(void)
{
        return default_create_format;
}
#endif

/* static */
NC *
ncmpii_dup_NC(const NC *ref)
{
        NC *ncp;

        ncp = (NC *) NCI_Malloc(sizeof(NC));
        if(ncp == NULL)
                return NULL;
        bzero(ncp, sizeof(NC));

        if(ncmpii_dup_NC_dimarrayV(&ncp->dims, &ref->dims) != NC_NOERR)
                goto err;
        if(ncmpii_dup_NC_attrarrayV(&ncp->attrs, &ref->attrs) != NC_NOERR)
                goto err;
        if(ncmpii_dup_NC_vararrayV(&ncp->vars, &ref->vars) != NC_NOERR)
                goto err;

        ncp->xsz = ref->xsz;
        ncp->begin_var = ref->begin_var;
        ncp->begin_rec = ref->begin_rec;
        ncp->recsize = ref->recsize;
        NC_set_numrecs(ncp, NC_get_numrecs(ref));
        return ncp;
err:
        ncmpii_free_NC(ncp);
        return NULL;
}


/*
 *  Verify that this is a user nc_type
 * Formerly
NCcktype()
 * Sense of the return is changed.
 */
int
ncmpii_cktype(nc_type type)
{
        switch((int)type){
        case NC_BYTE:
        case NC_CHAR:
        case NC_SHORT:
        case NC_INT:
        case NC_FLOAT:
        case NC_DOUBLE:
                return(NC_NOERR);
        }
        return(NC_EBADTYPE);
}


/*
 * How many objects of 'type'
 * will fit into xbufsize?
 */
MPI_Offset
ncmpix_howmany(nc_type type, MPI_Offset xbufsize)
{
        switch(type){
        case NC_BYTE:
        case NC_CHAR:
                return xbufsize;
        case NC_SHORT:
                return xbufsize/X_SIZEOF_SHORT;
        case NC_INT:
                return xbufsize/X_SIZEOF_INT;
        case NC_FLOAT:
                return xbufsize/X_SIZEOF_FLOAT;
        case NC_DOUBLE:
                return xbufsize/X_SIZEOF_DOUBLE;
        default:
                assert("ncmpix_howmany: Bad type" == 0);
                return(0);
        }
}

#define        D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

/*
 * Compute each variable's 'begin' offset,
 * update 'begin_rec' as well.
 */
static int
NC_begins(NC *ncp,
          MPI_Offset h_minfree, MPI_Offset v_align,
          MPI_Offset v_minfree, MPI_Offset r_align)
{
    MPI_Offset ii, sizeof_off_t, max_xsz;
    off_t index = 0;
    NC_var **vpp;
    NC_var *last = NULL;

    if (v_align == NC_ALIGN_CHUNK)  /* for non-record variables */
        v_align = ncp->chunk;
    if (r_align == NC_ALIGN_CHUNK)  /* for record variables */
        r_align = ncp->chunk;

    if (fIsSet(ncp->flags, NC_64BIT_OFFSET) ||
        fIsSet(ncp->flags, NC_64BIT_DATA)) {
        sizeof_off_t = 8;
    } else {
        sizeof_off_t = 4;
    }

    /* get the header size (un-aligned one) */
    ncp->xsz = ncmpii_hdr_len_NC(ncp, sizeof_off_t);

    if (ncp->vars.nelems == 0) /* no variable has been defined */
        return NC_NOERR;

    /* Since it is allowable to have incocnsistent header's metadata,
       header sizes, ncp->xsz, may be different among processes.
       Hence, we need to use the max size among all processes to make
       sure everyboday has the same size of reserved space for header */
    MPI_Allreduce(&ncp->xsz, &max_xsz, 1, MPI_LONG_LONG_INT, MPI_MAX,
                  ncp->nciop->comm);

    /* only (re)calculate begin_var if there is not sufficient space in header
       or start of non-record variables is not aligned as requested by valign */
    if (ncp->begin_var < max_xsz + h_minfree ||
        ncp->begin_var != D_RNDUP(ncp->begin_var, v_align) ) {
        index = (off_t) max_xsz;
        ncp->begin_var = D_RNDUP(index, v_align);
        if (ncp->begin_var < index + h_minfree)
            ncp->begin_var = D_RNDUP(index + (off_t)h_minfree, v_align);
    }

    index = ncp->begin_var;
    /* this is the aligned strating file offset for the first variable,
       also the header reserved size */

    /* Now find the starting file offsets for all variables.
       loop thru vars, first pass is for the 'non-record' vars */
    vpp = ncp->vars.value;
    for (ii=0; ii<ncp->vars.nelems ; ii++, vpp++) {
        if (IS_RECVAR(*vpp))
            /* skip record variables on this pass */
            continue;
#if 0
fprintf(stderr, "    VAR %lld %s: %lld\n", ii, (*vpp)->name->cp, index);
#endif
        if (sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0))
            return NC_EVARSIZE;

        /* this will pad out non-record variables with zero to the
         * requested alignment.  record variables are a bit trickier.
         * we don't do anything special with them */
        (*vpp)->begin = D_RNDUP(index, v_align);
        index = (*vpp)->begin + (*vpp)->len;
    }

    /* only (re)calculate begin_rec if there is not sufficient
       space at end of non-record variables or if start of record
       variables is not aligned as requested by r_align */
    if (ncp->begin_rec < index + v_minfree ||
        ncp->begin_rec != D_RNDUP(ncp->begin_rec, r_align) ) {
        ncp->begin_rec = D_RNDUP(index, r_align);
        if (ncp->begin_rec < index + v_minfree)
            ncp->begin_rec = D_RNDUP(index + (off_t)v_minfree, r_align);
    }
    index = ncp->begin_rec;

    ncp->recsize = 0;

    /* loop thru vars, second pass is for the 'record' vars */
    vpp = (NC_var **)ncp->vars.value;
    for (ii=0; ii<ncp->vars.nelems; ii++, vpp++) {
        if (!IS_RECVAR(*vpp))
            /* skip non-record variables on this pass */
            continue;

#if 0
fprintf(stderr, "    REC %lld %s: %lld\n", ii, (*vpp)->name->cp, index);
#endif
        if (sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0))
            return NC_EVARSIZE;

        /* A few attempts at aligning record variables have failed
         * (either with range error or 'value read not that expected',
         * or with an error in ncmpi_redef )).  Not sufficent to align
         * 'begin', but haven't figured out what else to adjust */
        (*vpp)->begin = index;
        index += (*vpp)->len;
        /* check if record size must fit in 32-bits */
#if SIZEOF_OFF_T == SIZEOF_SIZE_T && SIZEOF_SIZE_T == 4
        if (ncp->recsize > X_UINT_MAX - (*vpp)->len)
            return NC_EVARSIZE;
#endif
        ncp->recsize += (*vpp)->len;
        last = (*vpp);
    }

    /*
     * for special case of exactly one record variable, pack value
     */
    /* if there is exactly one record variable, then there is no need to
     * pad for alignment -- there's nothing after it */
    if (last != NULL && ncp->recsize == last->len)
        ncp->recsize = *last->dsizes * last->xsz;

    if (NC_IsNew(ncp))
        NC_set_numrecs(ncp, 0);

    return NC_NOERR;
}

#define NC_NUMRECS_OFFSET 4
#define NC_NUMRECS_EXTENT 4

/*
 * Read just the numrecs member.
 * (A relatively expensive way to do things.)
 */

 
int
ncmpii_read_numrecs(NC *ncp) {
  int status = NC_NOERR, mpireturn;
  MPI_Offset nrecs;
  void *buf, *pos;
  MPI_Status mpistatus;
  int rank;
  int sizeof_t;

  MPI_Comm_rank(ncp->nciop->comm, &rank);
 
  assert(!NC_indef(ncp));

  if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
                sizeof_t = X_SIZEOF_LONG;
  } else {
                 sizeof_t = X_SIZEOF_SIZE_T;
  }
 
  pos = buf = (void *)NCI_Malloc(sizeof_t);

  /* fileview is already entire file visible and pointer to zero */

  if (fIsSet(ncp->flags, NC_64BIT_DATA))
      mpireturn = MPI_File_read_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET+4,
                                   buf, sizeof_t, MPI_BYTE, &mpistatus);
  else
      mpireturn = MPI_File_read_at(ncp->nciop->collective_fh, NC_NUMRECS_OFFSET,
                                   buf, sizeof_t, MPI_BYTE, &mpistatus);
 
  if (mpireturn != MPI_SUCCESS) {
    ncmpii_handle_error(rank, mpireturn, "MPI_File_read_at");
    return NC_EREAD;
  } 

  status = ncmpix_get_size_t((const void **)&pos, &nrecs, sizeof_t);
  ncp->numrecs = nrecs;
 
  NCI_Free(buf);
 
  return status;
}
 
/*
 * Write out just the numrecs member.
 * (A relatively expensive way to do things.)
 */

/*
 * Collective operation implicit
 */

int
ncmpii_write_numrecs(NC *ncp)
{
    /* this function is only called by 3 places
       1) collective put APIs in mpinetcdf.c that write record variables
          and numrecs has change
       2) ncmpi_begin_indep_data() and ncmpi_end_indep_data()
       3) ncmpii_NC_sync() below
     */
    int rank, len, status = NC_NOERR, mpireturn;
    void *buf, *pos; 
    MPI_Offset nrecs;
    MPI_Status mpistatus;
    MPI_Comm comm;
  
    assert(!NC_readonly(ncp));
    assert(!NC_indef(ncp));
    /* this function is only called by APIs in data mode */

    comm = ncp->nciop->comm;
    MPI_Comm_rank(comm, &rank);

    /* there is a MPI_Allreduce() to ensure ncp->numrecs the same across all
       processes prior to this call */
    nrecs = ncp->numrecs;
    if (ncp->flags & NC_64BIT_DATA)
        len = X_SIZEOF_LONG;
    else
        len = X_SIZEOF_SIZE_T;
    pos = buf = (void *)NCI_Malloc(len);
    status = ncmpix_put_size_t(&pos, &nrecs, len);

    /* file view is already reset to entire file visible */
    if (rank == 0) {
        mpireturn = MPI_File_write_at(ncp->nciop->collective_fh,
                                      NC_NUMRECS_OFFSET, buf, len,
                                      MPI_BYTE, &mpistatus); 
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_handle_error(rank, mpireturn, "MPI_File_write_at");
            return NC_EWRITE;
        }
    }

    fClr(ncp->flags, NC_NDIRTY);  

    NCI_Free(buf);

    return status;
}

/*
 * Read in the header
 * It is expensive.
 */

int
ncmpii_read_NC(NC *ncp) {
  int status = NC_NOERR;

  ncmpii_free_NC_dimarrayV(&ncp->dims);
  ncmpii_free_NC_attrarrayV(&ncp->attrs);
  ncmpii_free_NC_vararrayV(&ncp->vars); 

  status = ncmpii_hdr_get_NC(ncp);

  if(status == NC_NOERR)
    fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

  return status;
}

/*
 * Write out the header
 */

static int
write_NC(NC *ncp)
{
    int status = NC_NOERR, mpireturn, rank;
    void *buf;
    MPI_Offset hsz; /* header size with 0-padding if needed */
    MPI_Status mpistatus;
 
    assert(!NC_readonly(ncp));
 
    MPI_Comm_rank(ncp->nciop->comm, &rank);

    if (NC_dofill(ncp)) { /* hsz is the header size with zero padding */
        /* we don't need this, as these paddings will never be accessed */
        hsz = MIN(ncp->begin_var, ncp->begin_rec);
        hsz = MAX(hsz, ncp->xsz);
    }
    else
        hsz = ncp->xsz;

    buf = (void *)NCI_Malloc(hsz); /* header buffer for I/O */
    if (hsz > ncp->xsz)
        bzero((char*)buf+ncp->xsz, hsz - ncp->xsz);

    status = ncmpii_hdr_put_NC(ncp, buf); /* copy header to buffer */
    if (status != NC_NOERR) {
        NCI_Free(buf);
        return status;
    }

    /* check the header consistency across all processes */
    status = NC_check_header(ncp->nciop->comm, buf, ncp->xsz, ncp);
    if (status != NC_NOERR) {
        NCI_Free(buf);
        return status;
    }

    /* the fileview is already entire file visible */

    /* only rank 0's header gets written to the file */
    if (rank == 0) {
        mpireturn = MPI_File_write_at(ncp->nciop->collective_fh, 0, buf,
                                      hsz, MPI_BYTE, &mpistatus);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_handle_error(rank, mpireturn, "MPI_File_write_at");
            return NC_EWRITE;
        }
    }
    fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);
    NCI_Free(buf);
 
    return status;
} 

/*
 * Write the header or the numrecs if necessary.
 */
int
ncmpii_NC_sync(NC  *ncp,
               int  doFsync)
{
    /* this function is called from four places:
       1) changing header by put APIs in data mode and NC_doHsync(ncp) is true
          and these put APIs are called in data mode
       2) ncmpi_sync()
       3) ncmpi_abort()
       4) ncmpii_NC_close()
       only 1) needs to call MPI_File_sync() here
     */
    int status = NC_NOERR, didWrite = 0;
    MPI_Offset mynumrecs, numrecs;

    assert(!NC_readonly(ncp));

    /* collect and set the max numrecs due to difference by independent write */
    mynumrecs = ncp->numrecs;
    MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX,
                  ncp->nciop->comm);
    if (numrecs > ncp->numrecs) {
        ncp->numrecs = numrecs;
        set_NC_ndirty(ncp); /* set numrecs is dirty flag */
    }

    if (NC_hdirty(ncp)) {  /* header is dirty */
        /* write_NC() will also write numrecs and clear NC_NDIRTY */
        status = write_NC(ncp);
        didWrite = 1;
    }
    else if (NC_ndirty(ncp)) {  /* numrecs is dirty */
        status = ncmpii_write_numrecs(ncp);
        didWrite = 1;
    }

    if (doFsync && didWrite)
        /* calling fsync only in data mode and header was just updataed above.
           If this API is called from ncmpi_sync() or ncmpi_abort(), fsync
           is no needed because a later fsync call is laready in those APIs */
        ncmpiio_sync(ncp->nciop);

    return status;
}


/*
 * header size increases, shift all record and non-record variables down
 */
static int
move_data_r(NC *ncp, NC *old) {
    /* no new record or non-record variable inserted, header size increases,
     * must shift (move) the whole contiguous data part down
     * Note ncp->numrecs may be > old->numrecs
     */
    return ncmpiio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
                        old->begin_rec - old->begin_var +
                        ncp->recsize * ncp->numrecs);
}

/*
 * Move the record variables down,
 * re-arrange records as needed
 * Fill as needed.
 */
static int
move_recs_r(NC *ncp, NC *old) {
    int status;
    MPI_Offset recno;
    const MPI_Offset nrecs = ncp->numrecs;
    const MPI_Offset ncp_recsize = ncp->recsize;
    const MPI_Offset old_recsize = old->recsize;
    const off_t ncp_off = ncp->begin_rec;
    const off_t old_off = old->begin_rec;

    assert(ncp_recsize >= old_recsize);

    if (ncp_recsize == old_recsize) {
        if (ncp_recsize == 0) /* no record variable defined yet */
            return NC_NOERR;

        /* No new record variable inserted, move all record variables as a whole */
        status = ncmpiio_move(ncp->nciop, ncp_off, old_off, ncp_recsize * nrecs);
        if (status != NC_NOERR)
            return status;
    } else {
        /* new record variables inserted, move one whole record at a time */
        for (recno = nrecs-1; recno >= 0; recno--) {
            status = ncmpiio_move(ncp->nciop, 
                                  ncp_off+recno*ncp_recsize, 
                                  old_off+recno*old_recsize, 
                                  old_recsize);
            if (status != NC_NOERR)
                return status;
        } 
    }
  
    return NC_NOERR;
}


/*
 * Move the "non record" variables "out". 
 * Fill as needed.
 */

static int
move_vars_r(NC *ncp, NC *old) {
  return ncmpiio_move(ncp->nciop, ncp->begin_var, old->begin_var, 
                   old->begin_rec - old->begin_var); 
}

#if 0
 /*
 * Given a valid ncp, return NC_EVARSIZE if any variable has a bad len 
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
static int
ncmpii_NC_check_vlens(NC *ncp)
{
    NC_var **vpp;
    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes.  This is different for format 1
       and format 2. */
    MPI_Offset vlen_max;
    MPI_Offset ii;
    MPI_Offset large_vars_count;
    MPI_Offset rec_vars_count;
    int last=-1;

    if(ncp->vars.nelems == 0) 
       return NC_NOERR;

    if ((ncp->flags & NC_64BIT_OFFSET) && sizeof(off_t) > 4) {
       /* CDF2 format and LFS */
       vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    } else if ((ncp->flags & NC_64BIT_DATA) && sizeof(off_t) > 4) {
       vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    } else {
       /* CDF1 format */
       vlen_max = X_INT_MAX - 3;
    }
    /* Loop through vars, first pass is for non-record variables.   */
    large_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
       if( !IS_RECVAR(*vpp) ) {
           last = 0;
           if( ncmpii_NC_check_vlen(*vpp, vlen_max) == 0 ) {
               large_vars_count++;
               last = 1;
           }
       } else {
         rec_vars_count++;
       }
    }
    /* OK if last non-record variable size too large, since not used to 
       compute an offset */
    if( large_vars_count > 1) { /* only one "too-large" variable allowed */
      return NC_EVARSIZE;
    }
    /* and it has to be the last one */ 
    if( large_vars_count == 1 && last == 0) { 
      return NC_EVARSIZE;
    }
    if( rec_vars_count > 0 ) {
       /* and if it's the last one, there can't be any record variables */
       if( large_vars_count == 1 && last == 1) {
           return NC_EVARSIZE;
       }
       /* Loop through vars, second pass is for record variables.   */
       large_vars_count = 0;
       vpp = ncp->vars.value;
       for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
           if( IS_RECVAR(*vpp) ) {
               last = 0;
               if( ncmpii_NC_check_vlen(*vpp, vlen_max) == 0 ) {
                   large_vars_count++;
                   last = 1;
               }
           }
       }
       /* OK if last record variable size too large, since not used to 
          compute an offset */
       if( large_vars_count > 1) { /* only one "too-large" variable allowed */
           return NC_EVARSIZE;
       }
       /* and it has to be the last one */ 
       if( large_vars_count == 1 && last == 0) { 
           return NC_EVARSIZE;
       }
    }
    return NC_NOERR;
}
#endif
 
int 
ncmpii_NC_enddef(NC *ncp) {
    int flag=0, alignment=512, status = NC_NOERR;
    char value[MPI_MAX_INFO_VAL];

    assert(!NC_readonly(ncp));
    /* assert(NC_indef(ncp));   This has alredy been checked */

    /* TODO: create a new member in the header struct for hint
             nc_file_header_alignment. This hint is passed from user's
             MPI_Info object at open time and will only be used by pnetcdf
             library. Note if this hint is passed to MPI-IO, it will be
             lost. So, it must be retrieved during file open.
             For example, ncp->nciop->info->header_alignment

       Currently, we use striping_unit or 512 for it.
     */
    ncmpiio_get_hint(ncp, "striping_unit", value, &flag);
  
    if (flag) 
        alignment=atoi(value);
    /* negative or zero alignment? can't imagine what that would even mean.
     * retain default 512 byte alignment (as in serial netcdf) */
    if (alignment <= 0)
        alignment = 512;

    /* NC_begins: pnetcdf doesn't expose an equivalent to nc__enddef, but we can
     * acomplish the same thing with calls to NC_begins
     * to compute each variable's 'begin' offset */
    NC_begins(ncp, 0, alignment, 0, alignment);
 
    /* serial netcdf calls a check on dimension lenghths here */
    /* To be updated */

    if (ncp->old != NULL) {
        /* the current define mode was entered from ncmpi_redef,
           not from ncmpi_create */

        assert(!NC_IsNew(ncp));
        assert(fIsSet(ncp->flags, NC_INDEF));
        assert(ncp->begin_rec >= ncp->old->begin_rec);
        assert(ncp->begin_var >= ncp->old->begin_var);
        assert(ncp->vars.nelems >= ncp->old->vars.nelems);
        /* ncp->numrecs has already sync-ed in ncmpi_redef */

        if (ncp->vars.nelems != 0) { /* number of record and non-record variables */
            if (ncp->begin_var > ncp->old->begin_var) {
                /* header size increases, shift the entire data part down */
                /* shift record variables first */
                status = move_recs_r(ncp, ncp->old);
                if (status != NC_NOERR)
                    return status;
                /* shift non-record variables */
                status = move_vars_r(ncp, ncp->old);
                if (status != NC_NOERR)
                    return status;
            }
            else if (ncp->begin_rec > ncp->old->begin_rec ||
                     ncp->recsize   > ncp->old->recsize) {
                /* number of non-record variables increases, or
                   number of records of record variables increases,
                   shift and move all record variables down */
                status = move_recs_r(ncp, ncp->old);
                if (status != NC_NOERR)
                    return status;
           }
       }
    } /* ... ncp->old != NULL */
 
    /* write header to file, also sync header buffer across all processes */
    status = write_NC(ncp);
    if (status != NC_NOERR)
        return status;
 
    if (ncp->old != NULL) {
        ncmpii_free_NC(ncp->old);
        ncp->old = NULL;
    }
    fClr(ncp->flags, NC_CREAT | NC_INDEF);

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE))
        /* calling MPI_File_sync() */
        ncmpiio_sync(ncp->nciop);

    return NC_NOERR;
}

#if 0
static int 
enddef(NC *ncp)
{
  assert(!NC_readonly(ncp));
  if(!NC_indef(ncp))
    return(NC_ENOTINDEFINE);

  NC_begins(ncp, 0, 1, 0, 1);

  if(ncp->old != NULL)
  {
    ncmpii_free_NC(ncp->old);
    ncp->old = NULL;
  }

  fClr(ncp->flags, NC_CREAT | NC_INDEF);

  return NC_NOERR;
}
#endif


/* Public */

int 
ncmpii_NC_close(NC *ncp) {
    int num_reqs, status=NC_NOERR, *req_ids=NULL, *statuses=NULL;
    NC_req *cur_req;

    if (NC_indef(ncp)) { /* currently in define mode */
        status = ncmpii_NC_enddef(ncp); /* TODO: defaults */
        if (status != NC_NOERR ) {
            /* To do: Abort new definition, if any */
            if (ncp->old != NULL) {
                ncmpii_free_NC(ncp->old);
                ncp->old = NULL;
                fClr(ncp->flags, NC_INDEF);
            }
        }
    }

    /* cancel or complete all outstanding nonblocking I/O */
    num_reqs = 0;
    cur_req = ncp->head;
    while (cur_req != NULL) {
        num_reqs++;
        cur_req = cur_req->next;
    }
    if (num_reqs > 0) { /* fill in req_ids[] */
        req_ids = (int*) NCI_Malloc(2 * num_reqs * sizeof(int));
        statuses = req_ids + num_reqs;
        num_reqs = 0;
        cur_req = ncp->head;
        while (cur_req != NULL) {
            req_ids[num_reqs++] = cur_req->id;
            cur_req = cur_req->next;
        }
    }
#if COMPLETE_NONBLOCKING_IO
    ncmpii_wait(ncp, COLL_IO, num_reqs, req_ids, statuses);
#else
    if (num_reqs > 0) ncmpii_cancel(ncp, num_reqs, req_ids, statuses);
#endif
    if (num_reqs > 0)
        NCI_Free(req_ids);

    if (!NC_readonly(ncp)) { /* file is open for write */
        /* check if header is dirty, if yes, flush it to file */
        status = ncmpii_NC_sync(ncp, 0);
        if (status != NC_NOERR)
            return status;
    }
 
    if (fIsSet(ncp->nciop->ioflags, NC_SHARE))
        /* calling MPI_File_sync() */
        ncmpiio_sync(ncp->nciop);

    ncmpiio_close(ncp->nciop, 0);
    ncp->nciop = NULL;
 
    ncmpii_del_from_NCList(ncp);
 
    ncmpii_free_NC(ncp);
 
    return status;
}


int
ncmpi_inq(int  ncid,
          int *ndimsp,
          int *nvarsp,
          int *nattsp,
          int *xtendimp)
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR)
        return status;

    if (ndimsp != NULL)
        *ndimsp = (int) ncp->dims.nelems;
    if (nvarsp != NULL)
        *nvarsp = (int) ncp->vars.nelems;
    if (nattsp != NULL)
        *nattsp = (int) ncp->attrs.nelems;
    if (xtendimp != NULL)
        *xtendimp = ncmpii_find_NC_Udim(&ncp->dims, NULL);

    return NC_NOERR;
}

int
ncmpi_inq_version(int ncid, int *NC_mode)
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if(status != NC_NOERR)
        return status;
        

    if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
        *NC_mode = NC_64BIT_DATA;
    } else if (fIsSet(ncp->flags, NC_64BIT_OFFSET)) {
        *NC_mode = NC_64BIT_OFFSET;
    } else {
        *NC_mode = 0;
    }
    return 0;
}



int 
ncmpi_inq_ndims(int ncid, int *ndimsp)
{
        int status;
        NC *ncp;

        status = ncmpii_NC_check_id(ncid, &ncp); 
        if(status != NC_NOERR)
                return status;

        if(ndimsp != NULL)
                *ndimsp = (int) ncp->dims.nelems;

        return NC_NOERR;
}

int 
ncmpi_inq_nvars(int ncid, int *nvarsp)
{
        int status;
        NC *ncp;

        status = ncmpii_NC_check_id(ncid, &ncp); 
        if(status != NC_NOERR)
                return status;

        if(nvarsp != NULL)
                *nvarsp = (int) ncp->vars.nelems;

        return NC_NOERR;
}

int 
ncmpi_inq_natts(int ncid, int *nattsp)
{
        int status;
        NC *ncp;

        status = ncmpii_NC_check_id(ncid, &ncp); 
        if(status != NC_NOERR)
                return status;

        if(nattsp != NULL)
                *nattsp = (int) ncp->attrs.nelems;

        return NC_NOERR;
}

int 
ncmpi_inq_unlimdim(int ncid, int *xtendimp)
{
        int status;
        NC *ncp;

        status = ncmpii_NC_check_id(ncid, &ncp); 
        if(status != NC_NOERR)
                return status;

        if(xtendimp != NULL)
                *xtendimp = ncmpii_find_NC_Udim(&ncp->dims, NULL);

        return NC_NOERR;
}
/*ARGSUSED*/

