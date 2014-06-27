/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include <mpi.h>

#include "nc.h"
#include "ncx.h"
#include "macro.h"
#ifdef ENABLE_SUBFILING
#include "subfile.h"
#endif

/* Prototypes for functions used only in this file */
static int ncmpii_end_indep_data(NC *ncp);


/* The const string below is for the RCS ident(1) command to find a string like
 * "\044Id: \100(#) PnetCDF library version 1.4.0 of 16 Nov 2013 $"
 * in the library file (libpnetcdf.a).
 */
static const char pnetcdf_libvers[] =
        "\044Id: \100(#) PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE" $";

/*----< ncmpi_inq_libvers() >------------------------------------------------*/
inline const char*
ncmpi_inq_libvers(void) {

    /* match the style used by netCDF API nc_inq_libvers()
     * for example, "4.3.0 of Jun 16 2013 12:11:30 $" */
    return PNETCDF_VERSION " of " PNETCDF_RELEASE_DATE;
}

/* Begin Of Dataset Functions */

/*----< ncmpi_create() >-----------------------------------------------------*/
/**  \ingroup datasets
Create a new netCDF file.

This function creates a new netCDF dataset, returning a netCDF ID that can
subsequently be used to refer to the netCDF dataset in other PnetCDF function
calls. The new netCDF dataset opened for write access and placed in define
mode, ready for you to add dimensions, variables, and attributes.

\param comm The MPI communicator. This API is a collective routine: all
processes must provide the same value for cmode, and all processes must provide
filenames that reference the same file. (Values for info may vary.) comm must
be an MPI intracommunicator. 

\param path The file name of the new netCDF dataset.

\param cmode The creation mode flag. The following flags are available:
NC_NOCLOBBER, NC_SHARE, NC_64BIT_OFFSET, and NC_64BIT_DATA.

\param info MPI info object. It is used to provide file access hints,including
existing MPI hints as well as PnetCDF hints.  For MPI hints, users are referred
to MPI user guide for further information. For PnetCDF hints see below.

\param ncidp Pointer to location where returned netCDF ID is to be stored.

<h2>The cmode Flag</h2>

The cmode flag is used to control the type of file created, and some aspects of
how it may be used. 

Setting NC_NOCLOBBER means you do not want to clobber (overwrite) an existing
dataset; an error (NC_EEXIST) is returned if the specified dataset already
exists.

The NC_SHARE flag is appropriate to synchronize the metadata across all
processes.  The metadata is store in the file header. For example, the number
of records may be increased by requests from only a subset of processes. Using
this flag will ensure the number of records (and other metadata) synchronized
across all the processes opening the same shared file.

Setting NC_64BIT_OFFSET causes PnetCDF to create a 64-bit offset format file
(CDF-2), instead of a netCDF classic format file.  The 64-bit offset format
imposes far fewer restrictions on large (i.e. over 2 GB) data files.  See Large
File Support (The PnetCDF Users Guide).

Setting NC_64BIT_DATA causes PnetCDF to create a 64-bit data format file
(CDF-5). The 64-bit data format allows define variables with more than 4
billion array elements. See Large File Support (The PnetCDF Users Guide).

A zero value (defined for convenience as NC_CLOBBER) specifies the default
behavior: overwrite any existing dataset with the same file name.

<h2>The info object Flag</h2>

Starting from version 1.3.1, the following PnetCDF hints are available:

- nc_header_align_size: This hint allows some extra space between the end of
the header describing the entire file and the first variable. If you have an
application that periodically wishes to add more variables to an already
existing file, expanding the file header size may result in an expensive move
of the entire data file to make room for the definition of the new variables.
Hence, setting this hint to a value that is big enough to accommodate any
additional variables means you may leave your application code as-is and yet
still see tremendous performance improvements.

- nc_var_align_size: If you are writing to a block-based parallel file system,
such as IBM's GPFS or Lustre, then an application write becomes a block write
at the file system layer. If a write straddles two blocks, then locks must be
acquired for both blocks. Aligning the start of a variable to a block boundary,
combined with collective I/O optimizations in the MPI-IO library can often
eliminate all unaligned file system accesses.

- nc_record_align_size: This hint aligns the starting file offset of the
record variable section.

- nc_header_read_chunk_size: PnetCDF reads the file headers in chunks. This
hint indicates the chunk size (in bytes). The default is 256 KB.

\returns ::NC_NOERR No error.
\returns ::NC_ENOMEM System out of memory.
\returns ::NC_EEXIST Specified file name exists when using NC_NOCLOBBER.
Can be use to check if the file exists.
\returns ::NC_EMULTIDEFINE_OMODE Bad file create/open mode or modes are
inconsistent across processes
\returns ::NC_EOFILE: Can not open/create file (MPI-IO errors)
\returns ::NC_EFILE: Unknown error in file operation

<h1>Examples</h1>

In this example we create a netCDF dataset named foo.nc; we want the dataset to
be created in the current directory only if a dataset with that name does not
already exist:

@code
     #include <mpi.h>
     #include <pnetcdf.h>
        ...
     int status;
     int ncid;
     MPI_Info info;
        ...
     MPI_Info_create (&info);
     MPI_Info_set (info, "romio_no_indep_rw",    "true");
     MPI_Info_set (info, "nc_header_align_size", "4194304");
     MPI_Info_set (info, "nc_var_align_size",    "1048576");
     MPI_Info_set (info, "nc_record_align_size", "1048576");

     status = ncmpi_create(MPI_COMM_WORLD, "foo.nc", NC_NOCLOBBER, info, &ncid);
     if (status != NC_NOERR) handle_error(status);
     MPI_Info_free(&info);
@endcode

In this example we create a netCDF dataset named foo.nc. It will
be in the CDF-5 format.

@code
     #include <mpi.h>
     #include <pnetcdf.h>
        ...
     int status;
     int ncid;
     int cmode = NC_NOCLOBBER | NC_64BIT_DATA;
     MPI_Info info = MPI_INFO_NULL;
        ...
     status = ncmpi_create(MPI_COMM_WORLD, "foo.nc", cmode, info, &ncid);
     if (status != NC_NOERR) handle_error(status);
@endcode
*/
int 
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    int err, status=NC_NOERR, safe_mode=0;
    char *env_str=NULL, *hint_str;
    MPI_Info   env_info;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp;

#ifdef PNC_DEBUG
    safe_mode = 1;
    /* this configure time setting will be overwritten by the run-time
     * environment variable PNETCDF_SAFE_MODE */
#endif
    /* get environment variable PNETCDF_SAFE_MODE
     * if it is set to 1, then we perform a strict parameter consistent test
     */
    env_str = getenv("PNETCDF_SAFE_MODE");
    if (env_str != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
    }

    if (safe_mode) {
        /* check if cmode is consistent with root's */
        int root_cmode=cmode;

        MPI_Bcast(&root_cmode, 1, MPI_INT, 0, comm);

        if (root_cmode != cmode) {
            int rank;
            MPI_Comm_rank(comm, &rank);
            /* cmodes are inconsistent, overwrite local cmode with root's */
            printf("rank %d: Warning - inconsistent file create mode, overwrite with root's\n",rank);
            cmode = root_cmode;
            status = NC_EMULTIDEFINE_OMODE;
        }
        /* when safe_mode is disabled, NC_EMULTIDEFINE_OMODE will be reported at
         * the time ncmpi_enddef() returns */
    }

    /* take hints from the environment variable PNETCDF_HINTS
     * a string of hints separated by ";" and each hint is in the
     * form of hint=value. E.g. cb_nodes=16;cb_config_list=*:6
     * If this environment variable is set, it  overrides any values that
     * were set by using calls to MPI_Info_set in the application code.
     */
    env_str = getenv("PNETCDF_HINTS");
    env_info = info;
    if (env_str != NULL) {
        if (info == MPI_INFO_NULL) {
            err = MPI_Info_create(&env_info);
            if (err != MPI_SUCCESS) /* ignore this error */
                env_info = MPI_INFO_NULL;
        }
        hint_str = strtok(env_str, ";");
        while (hint_str != NULL && env_info != MPI_INFO_NULL) {
            char key[128], *val;
            strcpy(key, hint_str);
            val = strchr(key, '=');
            *val = '\0';
            val++;
            /* printf("env hint: key=%s val=%s\n",key,val); */
            MPI_Info_set(env_info, key, val); /* override */
            hint_str = strtok(NULL, ";");
        }
    }

    /* get header chunk size from user info */
    if (env_info != MPI_INFO_NULL) {
        char value[MPI_MAX_INFO_VAL];
        int  flag;
        MPI_Info_get(env_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) chunksize = atoll(value);
    }

    /* allocate buffer for header object NC */
    if ((ncp = ncmpii_new_NC(&chunksize)) == NULL) 
        return NC_ENOMEM;

    ncp->safe_mode = safe_mode;
    ncp->old       = NULL;
#ifdef ENABLE_SUBFILING
    ncp->ncid_sf = -1; /* subfile ncid; init to -1 */ 
    ncp->nc_num_subfiles = 0; /* num_subfiles; init to 0 */
#endif
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

    err = ncmpiio_create(comm, path, cmode, env_info, ncp);  
    if (err != NC_NOERR) {
        ncmpii_free_NC(ncp);
        return err;
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

    if (env_info != info) MPI_Info_free(&env_info);

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
    int err, status=NC_NOERR, safe_mode=0;
    char *env_str=NULL, *hint_str;
    MPI_Info   env_info;
    MPI_Offset chunksize=NC_DEFAULT_CHUNKSIZE;
    NC *ncp;
  
#ifdef PNC_DEBUG
    safe_mode = 1;
    /* this configure time setting will be overwritten by the run-time
     * environment variable PNETCDF_SAFE_MODE */
#endif
    /* get environment variable PNETCDF_SAFE_MODE
     * if it is set to 1, then we perform a strict parameter consistent test
     */
    env_str = getenv("PNETCDF_SAFE_MODE");
    if (env_str != NULL) {
        if (*env_str == '0') safe_mode = 0;
        else                 safe_mode = 1;
    }

    if (safe_mode) {
        /* check if omode is consistent with root's */
        int root_omode=omode;

        /* Note if omode contains NC_NOWRITE, it is equivalent to NC_CLOBBER.
           In pnetcdf.h, they both are defined the same value, 0.
         */

        MPI_Bcast(&root_omode, 1, MPI_INT, 0, comm);

        if (root_omode != omode) {
            int rank;
            MPI_Comm_rank(comm, &rank);
            /* omodes are inconsistent, overwrite local omode with root's */
            printf("rank %d: Warning - inconsistent file open mode, overwrite with root's\n",rank);
            omode = root_omode;
            status = NC_EMULTIDEFINE_OMODE;
        }
    }

    /* take hints from the environment variable PNETCDF_HINTS
     * a string of hints separated by ";" and each hint is in the
     * form of hint=value. E.g. cb_nodes=16;cb_config_list=*:6
     * If this environment variable is set, it  overrides any values that
     * were set by using calls to MPI_Info_set in the application code.
     */
    env_str = getenv("PNETCDF_HINTS");
    env_info = info;
    if (env_str != NULL) {
        if (info == MPI_INFO_NULL) {
            err = MPI_Info_create(&env_info);
            if (err != MPI_SUCCESS) /* ignore this error */
                env_info = MPI_INFO_NULL;
        }
        hint_str = strtok(env_str, ";");
        while (hint_str != NULL && env_info != MPI_INFO_NULL) {
            char key[128], *val;
            strcpy(key, hint_str);
            val = strchr(key, '=');
            *val = '\0';
            val++;
            /* printf("env hint: key=%s val=%s\n",key,val); */
            MPI_Info_set(env_info, key, val); /* override */
            hint_str = strtok(NULL, ";");
        }
    }

    /* get header chunk size from user info, if provided */
    if (env_info != MPI_INFO_NULL) {
        char value[MPI_MAX_INFO_VAL];
        int  flag;
        MPI_Info_get(env_info, "nc_header_read_chunk_size", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        if (flag) chunksize = atoll(value);
    }

    ncp = ncmpii_new_NC(&chunksize);
    if (ncp == NULL)
        return NC_ENOMEM;

    ncp->safe_mode = safe_mode;
    ncp->old       = NULL;
#ifdef ENABLE_SUBFILING
    ncp->ncid_sf   = -1;
    ncp->nc_num_subfiles = 0;
#endif

    err = ncmpiio_open(comm, path, omode, env_info, ncp);
    if (err != NC_NOERR) {
        ncmpii_free_NC(ncp);
        return err;
    } 

    assert(ncp->flags == 0); 

    if (fIsSet(ncp->nciop->ioflags, NC_SHARE)) {
        /*
         * NC_SHARE implies sync up the number of records as well.
         * Note that other header changes are not shared
         * automatically.  Some sort of IPC (external to this package)
         * would be used to trigger a call to ncmpi_sync().
         */ 
        fSet(ncp->flags, NC_NSYNC);  /* sync numrecs */
        fSet(ncp->flags, NC_HSYNC);  /* sync header */
    }

    err = ncmpii_hdr_get_NC(ncp); /* read header from file */
    if (err != NC_NOERR) {
        ncmpiio_close(ncp->nciop, 0);
        ncmpii_free_NC(ncp);
        return err;
    }
    ncp->head = NULL;
    ncp->tail = NULL;

    ncmpii_add_to_NCList(ncp);
    *ncidp = ncp->nciop->fd;

#ifdef ENABLE_SUBFILING
    /* check attr for subfiles */
    nc_type type;
    MPI_Offset attlen;
    int ndims1, nvars1, natts1, unlimdimid1;
    int i;
    /* we report the first encountered status if there is an error */
    err = ncmpi_inq_att(ncp->nciop->fd, NC_GLOBAL, "num_subfiles",
                        &type, &attlen);
    if (err == NC_NOERR) {
        err = ncmpi_get_att_int(ncp->nciop->fd, NC_GLOBAL, "num_subfiles",
                                &ncp->nc_num_subfiles); 
        if (status == NC_NOERR) status = err;
        /* TODO: check err */

        err = ncmpi_inq(ncp->nciop->fd, &ndims1, &nvars1, &natts1, &unlimdimid1);
        if (status == NC_NOERR) status = err;

        for (i=0; i<nvars1; i++) {
            err = ncmpi_get_att_int(ncp->nciop->fd, i, "num_subfiles",
                                    &ncp->vars.value[i]->num_subfiles); 
            if (status == NC_NOERR) status = err;

            if (ncp->vars.value[i]->num_subfiles > 1) {
                err = ncmpi_get_att_int(ncp->nciop->fd, i, "ndims_org",
                                        &ncp->vars.value[i]->ndims_org);
                if (status == NC_NOERR) status = err;
            }
        }
    }

    if (ncp->nc_num_subfiles > 1) {
        err = ncmpii_subfile_open(ncp, &ncp->ncid_sf);
        if (status == NC_NOERR) status = err;
    }
#endif

    if (env_info != info) MPI_Info_free(&env_info);

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

/*----< ncmpi_inq_file_info() >-----------------------------------------------*/
int
ncmpi_inq_file_info(int       ncid,
                    MPI_Info *info_used)
{
    int mpireturn, status=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

#ifdef HAVE_MPI_INFO_DUP
    mpireturn = MPI_Info_dup(ncp->nciop->mpiinfo, info_used);
    if (mpireturn != MPI_SUCCESS) {
        ncmpii_handle_error(mpireturn, "MPI_Info_dup");
        return NC_EFILE;
    }
#else
    mpireturn = MPI_File_get_info(ncp->nciop->collective_fh, info_used);
    if (mpireturn != MPI_SUCCESS) {
        ncmpii_handle_error(mpireturn, "MPI_File_get_info");
        return NC_EFILE;
    }
#endif

    sprintf(value, "%lld", ncp->nciop->hints.h_align);
    MPI_Info_set(*info_used, "nc_header_align_size", value);

    sprintf(value, "%lld", ncp->nciop->hints.v_align);
    MPI_Info_set(*info_used, "nc_var_align_size", value);

    sprintf(value, "%lld", ncp->nciop->hints.r_align);
    MPI_Info_set(*info_used, "nc_record_align_size", value);

    sprintf(value, "%lld", ncp->nciop->hints.header_read_chunk_size);
    MPI_Info_set(*info_used, "nc_header_read_chunk_size", value);

#ifdef ENABLE_SUBFILING
    sprintf(value, "%d", ncp->nciop->hints.num_subfiles);
    MPI_Info_set(*info_used, "nc_num_subfiles", value);
#endif

    return NC_NOERR;
}

/*----< ncmpi_get_file_info() >-----------------------------------------------*/
int
ncmpi_get_file_info(int       ncid,
                    MPI_Info *info_used)
{
    return ncmpi_inq_file_info(ncid, info_used);
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
        /* if open mode is inconsistent, then this return might cause parallel
         * program to hang */

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
        MPI_Allreduce(&mynumrecs, &numrecs, 1, MPI_OFFSET, MPI_MAX,
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

/*----< ncmpii_sync_numrecs() >-----------------------------------------------*/
/* synchronize the number of records in memory and file header, if the number
 * is increased
 * This function is called by collective APIs only
 */
int
ncmpii_sync_numrecs(NC         *ncp,
                    MPI_Offset  new_numrecs)
{
    /* Only put APIs and record variables reach this function */
    MPI_Offset max_numrecs;

    /* sync numrecs in memory across all processes */
    MPI_Allreduce(&new_numrecs, &max_numrecs, 1, MPI_OFFSET, MPI_MAX,
                  ncp->nciop->comm);

    /* let root process write numrecs to file header
     * Note that we must update numrecs in file header whenever a
     * collective write call to a record variable is completed, instead
     * of waiting until file close. This is in case when the program
     * aborts before reaching ncmpi_close(), the numrecs value in the
     * file header can be up-to-date.
     * ncmpii_write_numrecs() also updates ncp->numrecs
     */
    return ncmpii_write_numrecs(ncp, max_numrecs, NC_ndirty(ncp));
}
 
/*----< ncmpi_begin_indep_data() >-------------------------------------------*/
int
ncmpi_begin_indep_data(int ncid)
{
#ifndef DISABLE_FILE_SYNC
    int mpireturn;
#endif
    int status=NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    if (NC_indef(ncp))  /* must not be in define mode */
        return NC_EINDEFINE;

    if (NC_indep(ncp))  /* already in indep data mode */
        return NC_NOERR;

#ifndef DISABLE_FILE_SYNC
    if (!NC_readonly(ncp) && NC_collectiveFhOpened(ncp->nciop)) {
        /* MPI_File_sync() is collective */
        mpireturn = MPI_File_sync(ncp->nciop->collective_fh);
        if (mpireturn != MPI_SUCCESS) {
            ncmpii_handle_error(mpireturn, "MPI_File_sync");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = NC_EFILE;
        }
    }
#endif

    fSet(ncp->flags, NC_INDEP);

    return status;
}

/*----< ncmpi_end_indep_data() >---------------------------------------------*/
int 
ncmpi_end_indep_data(int ncid) {
    int status;
    NC *ncp;
 
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR) return status;

    if (!NC_indep(ncp)) /* must be in independent data mode */
        return NC_ENOTINDEP;

    return ncmpii_end_indep_data(ncp);
}

/*----< ncmpii_end_indep_data() >--------------------------------------------*/
static int 
ncmpii_end_indep_data(NC *ncp) {
#ifndef DISABLE_FILE_SYNC
    int mpireturn;
#endif
    int status=NC_NOERR;

    if (!NC_readonly(ncp)) {
        /* do memory and file sync for numrecs, number or records */
        status = ncmpii_sync_numrecs(ncp, ncp->numrecs);

#ifndef DISABLE_FILE_SYNC
        /* calling file sync for those already open the file */
        if (NC_independentFhOpened(ncp->nciop)) {
            /* MPI_File_sync() is collective */
            mpireturn = MPI_File_sync(ncp->nciop->independent_fh);
            if (mpireturn != MPI_SUCCESS) {
                ncmpii_handle_error(mpireturn, "MPI_File_sync");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = NC_EFILE;
            }
        }
#endif
    }

    fClr(ncp->flags, NC_INDEP);

    return status;
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

    return ncmpii_enddef(ncp);
}

/*----< ncmpi__enddef() >-----------------------------------------------------*/
int
ncmpi__enddef(int        ncid,
              MPI_Offset h_minfree,
              MPI_Offset v_align,
              MPI_Offset v_minfree,
              MPI_Offset r_align)
{
    int status;
    NC *ncp;

    /* check if file ID ncid is valid */
    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR) return status;

    if (!NC_indef(ncp)) /* must currently in define mode */
        return NC_ENOTINDEFINE;

    return ncmpii__enddef(ncp, h_minfree, v_align, v_minfree, r_align);
}

/*----< ncmpi_sync_numrecs() >------------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note that numrecs is always sync-ed in collective mode
 */
int
ncmpi_sync_numrecs(int ncid) {
    int status = NC_NOERR;
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (NC_indef(ncp)) 
        return NC_EINDEFINE;

    /* syn numrecs in memory (and file if NC_SHARE is set) */
    return ncmpii_sync_numrecs(ncp, ncp->numrecs);
}

/*----< ncmpi_sync() >--------------------------------------------------------*/
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

    /* write header to file in ncmpii_NC_sync((), but don't call fsync now,
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
        if (mpireturn != MPI_SUCCESS) {
            ncmpiio_free(ncp->nciop);
            return ncmpii_handle_error(mpireturn, "MPI_File_open");
        }

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

/*----< ncmpi_inq_striping() >------------------------------------------------*/
/* return file (system) striping settings, striping size and count, if they are
 * available from MPI-IO hint. Otherwise, 0s are returned.
 */
int
ncmpi_inq_striping(int  ncid,
                   int *striping_size,
                   int *striping_count)
{
    int flag, status=NC_NOERR;
    char value[MPI_MAX_INFO_VAL];
    NC *ncp;

    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (striping_size != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_unit", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_size = 0;
        if (flag) *striping_size = atoi(value);
    }

    if (striping_count != NULL) {
        MPI_Info_get(ncp->nciop->mpiinfo, "striping_factor", MPI_MAX_INFO_VAL-1,
                     value, &flag);
        *striping_count = 0;
        if (flag) *striping_count = atoi(value);
    }
    return NC_NOERR;
}

/*----< ncmpi_inq_malloc_size() >--------------------------------------------*/
/* report the current aggregate size allocated by malloc, yet to be freed */
int ncmpi_inq_malloc_size(MPI_Offset *size)
{
#ifdef PNC_DEBUG
    ncmpii_inq_malloc_size(size);
    return NC_NOERR;
#else
    return NC_ENOTENABLED;
#endif
}

/*----< ncmpi_inq_malloc_max_size() >----------------------------------------*/
/* get the max watermark ever researched by malloc (aggregated amount) */
int ncmpi_inq_malloc_max_size(MPI_Offset *size)
{
#ifdef PNC_DEBUG
    ncmpii_inq_malloc_max_size(size);
    return NC_NOERR;
#else
    return NC_ENOTENABLED;
#endif
}

/*----< ncmpi_inq_malloc_list() >--------------------------------------------*/
/* walk the malloc tree and print yet-to-be-freed malloc residues */
int ncmpi_inq_malloc_list(void)
{
#ifdef PNC_DEBUG
    ncmpii_inq_malloc_list();
    return NC_NOERR;
#else
    return NC_ENOTENABLED;
#endif
}

/*----< ncmpi_inq_files_opened() >-------------------------------------------*/
int
ncmpi_inq_files_opened(int *num, int *ncids)
{
    return ncmpii_inq_files_opened(num, ncids);
}

