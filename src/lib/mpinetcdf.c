/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include "macro.h"

/* Prototypes for functions used only in this file */
static int ncmpii_begin_indep_data(NC *ncp);
static int ncmpii_end_indep_data(NC *ncp);

/*----< ncmpi_inq_libvers() >------------------------------------------------*/
inline const char*
ncmpi_inq_libvers(void) {
    return "version = " PNETCDF_VERSION " of 24 Sep 2012";
}

/* Begin Of Dataset Functions */

/*----< ncmpi_create() >-----------------------------------------------------*/
int 
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    int status;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp;

    /* check if cmode are consistent across all processes */
    int my_cmode, cmode_sum, nprocs;
    MPI_Comm_size(comm, &nprocs);

    /* Note if cmode contains NC_NOWRITE, it is equivalent to NC_CLOBBER.
       In pnetcdf.h, they both are defined the same value, 0.
     */

    my_cmode = 1;
    if (cmode & NC_64BIT_OFFSET)  my_cmode = 2;
    if (cmode & NC_64BIT_DATA)    my_cmode = 5;

    MPI_Allreduce(&my_cmode, &cmode_sum, 1, MPI_INT, MPI_SUM, comm);
    if (cmode_sum != my_cmode * nprocs) {
        // fprintf(stderr,"Error: create modes are inconsistent\n");
        *ncidp = -1;  /* cause NC_EBADID for any further opertion */
        return NC_ECMODE;
    }
    /* if cmodes are inconsistent, then it is a fatal error to continue */

    /* get header chunk size from user info */
    if (info != MPI_INFO_NULL) {
        char value[MPI_MAX_INFO_VAL];
        int  flag;
        MPI_Info_get(info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) chunksize = atoll(value);
    }

    /* allocate buffer for header object NC */
    if ((ncp = ncmpii_new_NC(&chunksize)) == NULL) 
        return NC_ENOMEM;

    ncp->old = NULL;
    assert(ncp->flags == 0);

    /* set the file format version beased on the create mode, cmode */
    if (fIsSet(cmode, NC_64BIT_OFFSET)) {
        /* unlike serial netcdf, we will not bother to support
         * NC_64BIT_OFFSET on systems with off_t smaller than 8 bytes.
         * serial netcdf has proven it's possible if datasets are small, but
         * that's a hassle we don't want to worry about */
        if (sizeof(off_t) != 8)
            return NC_ESMALL;
        fSet(ncp->flags, NC_64BIT_OFFSET);
    } else if (fIsSet(cmode, NC_64BIT_DATA)) {
        if (sizeof(MPI_Offset) <  8)
            return NC_ESMALL;
        fSet(ncp->flags, NC_64BIT_DATA);
    } else {
        fSet(ncp->flags, NC_32BIT);
    }

    /* find the true header size (not-yet aligned) */
    ncp->xsz = ncmpii_hdr_len_NC(ncp);

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

    /* the linked list storing the outstanding non-blocking requests */
    ncp->head = NULL;
    ncp->tail = NULL;

    /* add to the linked list of opened files */
    ncmpii_add_to_NCList(ncp);
    *ncidp = ncp->nciop->fd;

    return status;
}

/*----< ncmpi_open() >-------------------------------------------------------*/
int
ncmpi_open(MPI_Comm    comm,
           const char *path,
           int         omode,
           MPI_Info    info,
           int        *ncidp)
{
    int status = NC_NOERR;
    NC *ncp;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
  
    /* get header chunk size from user info, if provided */
    if (info != MPI_INFO_NULL) {
        char value[MPI_MAX_INFO_VAL];
        int  flag;
        MPI_Info_get(info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) chunksize = atoll(value);
    }

    ncp = ncmpii_new_NC(&chunksize);
    if (ncp == NULL)
        return NC_ENOMEM;

    ncp->old = NULL;

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
    ncp->head = NULL;
    ncp->tail = NULL;

    ncmpii_add_to_NCList(ncp);
    *ncidp = ncp->nciop->fd;

    return status;
}

/*----< ncmpi_inq_format() >-------------------------------------------------*/
int
ncmpi_inq_format(int  ncid,
                 int *formatp) /* out */
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
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

/*----< ncmpi_inq_file_format() >--------------------------------------------*/
int
ncmpi_inq_file_format(char *filename,
                      int  *formatp) /* out */
{
    int ncid, status;
    NC *ncp;
        
    /* open file for reading its header */
    status = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL,
                        &ncid);
    if (status == NC_ENOTNC)
        *formatp = NC_FORMAT_UNKNOWN;
    if (status != NC_NOERR)
        return status;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
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

/*----< ncmpi_get_file_info() >----------------------------------------------*/
int
ncmpi_get_file_info(int       ncid,
                    MPI_Info *info_used)
{
    int mpireturn, status=NC_NOERR, mpi_err=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

#ifdef HAVE_MPI_INFO_DUP
    mpireturn = MPI_Info_dup(ncp->nciop->mpiinfo, info_used);
    CHECK_MPI_ERROR(mpireturn, "MPI_Info_dup", NC_EFILE);
#else
    mpireturn = MPI_File_get_info(ncp->nciop->collective_fh, info_used);
    CHECK_MPI_ERROR(mpireturn, "MPI_File_get_info", NC_EFILE);
#endif

    sprintf(value, "%lld", ncp->nciop->hints.header_align_size);
    MPI_Info_set(*info_used, "nc_header_align_size", value);

    sprintf(value, "%lld", ncp->nciop->hints.var_align_size);
    MPI_Info_set(*info_used, "nc_var_align_size", value);

    sprintf(value, "%lld", ncp->nciop->hints.header_read_chunk_size);
    MPI_Info_set(*info_used, "nc_header_read_chunk_size", value);

    /* make NC error higher priority than MPI error */
    return (status != NC_NOERR) ? status : mpi_err;
}

/*----< ncmpi_redef() >------------------------------------------------------*/
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
        ncmpii_end_indep_data(ncp);

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
        MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_LONG_LONG_INT, MPI_MAX,
                      ncp->nciop->comm);
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

/*----< ncmpii_update_numrecs() >--------------------------------------------*/
int
ncmpii_update_numrecs(NC         *ncp,
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
        int localChange=0, doChange;

        if (NC_ndirty(ncp))
            localChange = 1;

        MPI_Allreduce(&localChange, &doChange, 1, MPI_INT, MPI_MAX,
                      ncp->nciop->comm);

        if (doChange) {
            /* all proc must agree on numrecs because this func is collective */
            MPI_Allreduce(&newnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT,
                          MPI_MAX, ncp->nciop->comm);
            status = ncmpii_write_numrecs(ncp);
            if (status != NC_NOERR)
                return status;

            /* fsync to disk */
            ncmpiio_sync(ncp->nciop);
        }
    }
    return status;
}
 
/*----< ncmpi_begin_indep_data() >-------------------------------------------*/
int
ncmpi_begin_indep_data(int ncid) {
    int status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    if (NC_indef(ncp))  /* must not be in define mode */
        return NC_EINDEFINE;

#if 0
    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_EINDEP; /* Should we skip this error? */
#endif
 
    return ncmpii_begin_indep_data(ncp);  
}

/*----< ncmpii_begin_indep_data() >------------------------------------------*/
static int
ncmpii_begin_indep_data(NC *ncp) {
    int mpireturn, status=NC_NOERR, mpi_err=NC_NOERR;

    if (!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
        /* do memory and file sync for numrecs, number or records */
        MPI_Offset oldnumrecs = ncp->numrecs;
        MPI_Allreduce(&oldnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT,
                      MPI_MAX, ncp->nciop->comm);
        status = ncmpii_write_numrecs(ncp);

        /* MPI_File_sync() is collective */
        mpireturn = MPI_File_sync(ncp->nciop->collective_fh);
        CHECK_MPI_ERROR(mpireturn, "MPI_File_sync", NC_EFILE);
    }

    fSet(ncp->flags, NC_INDEP);

    /* make NC error higher priority than MPI error */
    return (status != NC_NOERR) ? status : mpi_err;  
}

/*----< ncmpi_end_indep_data() >---------------------------------------------*/
int 
ncmpi_end_indep_data(int ncid) {
    int status = NC_NOERR;
    NC *ncp;
 
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    if (!NC_indep(ncp)) return NC_ENOTINDEP;

    return ncmpii_end_indep_data(ncp);
}

/*----< ncmpii_end_indep_data() >--------------------------------------------*/
static int 
ncmpii_end_indep_data(NC *ncp) {
    int mpireturn, status=NC_NOERR, mpi_err=NC_NOERR;

    if (!NC_readonly(ncp)) {
        /* do memory and file sync for numrecs, number or records */
        MPI_Offset oldnumrecs = ncp->numrecs;
        MPI_Allreduce(&oldnumrecs, &ncp->numrecs, 1, MPI_LONG_LONG_INT,
                      MPI_MAX, ncp->nciop->comm);
        status = ncmpii_write_numrecs(ncp);

        /* calling file sync for those already open the file */
        if (NC_independentFhOpened(ncp->nciop)) {
            /* MPI_File_sync() is collective */
            mpireturn = MPI_File_sync(ncp->nciop->independent_fh);
            CHECK_MPI_ERROR(mpireturn, "MPI_File_sync", NC_EFILE);
        }
    }

    fClr(ncp->flags, NC_INDEP);

    /* make NC error higher priority than MPI error */
    return (status != NC_NOERR) ? status : mpi_err;
}

/*----< ncmpi_enddef() >-----------------------------------------------------*/
int
ncmpi_enddef(int ncid) {
    int status;
    NC *ncp;

    /* check if file ID ncid is valid */
    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR) return status;

    if (!NC_indef(ncp)) /* must currently in define mode */
        return NC_ENOTINDEFINE;

    return ncmpii_NC_enddef(ncp);
}

/*----< ncmpi_sync() >-------------------------------------------------------*/
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

    /* write header to file in cmpii_NC_sync((), but don't call fsync now,
       fsync will be called below in ncmpiio_sync() */
    status = ncmpii_NC_sync(ncp, 0);
    if (status != NC_NOERR)
        return status;

    /* calling MPI_File_sync() */
    return ncmpiio_sync(ncp->nciop);
}

/*----< ncmpi_abort() >------------------------------------------------------*/
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
    if (status != NC_NOERR)
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

/*----< ncmpi_close() >------------------------------------------------------*/
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

/*----< ncmpi_delete() >-----------------------------------------------------*/
/* ncmpi_delete:
 * doesn't do anything to release resources, so call ncmpi_close before calling
 * this function.
 *
 * filename: the name of the
 * file we will remove.  info: mpi info, in case underlying file system needs
 * hints.
 */
int
ncmpi_delete(char     *filename,
             MPI_Info  info)
{
    int status = NC_NOERR;
    status = MPI_File_delete(filename, info);
    if (status != MPI_SUCCESS)
        return NC_EFILE;
    return NC_NOERR;
}

/*----< ncmpi_set_fill() >---------------------------------------------------*/
/* ncmpi_set_fill:
 * not actually implemented.  Anything other than NC_NOFILL is not supported.
 * Many codes use NC_NOFILL anyway, so this just gets us more source-portable
 * with existings serial netcdf codes.   Also provides a placeholder if someday
 * someone wants to implement all of set_fill 
 */
int
ncmpi_set_fill(int  ncid,
               int  fillmode,
               int *old_mode_ptr)
{
    int status = NC_NOERR;
    if (fillmode != NC_NOFILL)
        status = NC_EINVAL;
    return status;
}
                
/* End Of Dataset Functions */

/*----< ncmpii_check_mpifh() >-----------------------------------------------*/
int
ncmpii_check_mpifh(NC       *ncp,
                   MPI_File *mpifh,
                   MPI_Comm  comm,
                   int       collective)
{
    if (collective && NC_indep(ncp)) /* collective handle but in indep mode */
        return NC_EINDEP;

    if (!collective && !NC_indep(ncp)) /* indep handle but in collective mode */
        return NC_ENOTINDEP;

    if ( (collective && !NC_collectiveFhOpened(ncp->nciop))  ||
         (!collective && !NC_independentFhOpened(ncp->nciop)) ) {
  
        int mpireturn;
        mpireturn = MPI_File_open(comm, (char *)ncp->nciop->path,
                                  ncp->nciop->mpiomode, ncp->nciop->mpiinfo,
                                  mpifh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_check_mpi_file_open_error(ncp->nciop, mpireturn);

        if (collective)
            set_NC_collectiveFh(ncp->nciop);
        else
            set_NC_independentFh(ncp->nciop);
    }
    return NC_NOERR;
}

/*----< ncmpi_inq_put_size() >------------------------------------------------*/
/* returns the amount of writes, in bytes, committed to file system so far */
int 
ncmpi_inq_put_size(int         ncid,
                   MPI_Offset *size)
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    *size = ncp->nciop->put_size;
    return NC_NOERR;
}

/*----< ncmpi_inq_get_size() >------------------------------------------------*/
/* returns the amount of reads, in bytes, obtained from file system so far */
int 
ncmpi_inq_get_size(int         ncid,
                   MPI_Offset *size)
{
    int status;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    *size = ncp->nciop->get_size;
    return NC_NOERR;
}
