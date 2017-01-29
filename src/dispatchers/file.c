
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>   /* open() */
#include <unistd.h>  /* read(), close() */
#include <errno.h>   /* errno */

#include <dispatch.h>
#include <pnetcdf.h>

/* static variables are initialized to NULLs */
static PNC *pnc_filelist[NC_MAX_FILES];
static int  pnc_numfiles;

static int
add_to_PNCList(PNC *pncp,
               int *new_id)
{   
    int i;

    *new_id = -1;
    for (i=0; i<NC_MAX_FILES; i++) {
        if (pnc_filelist[i] == NULL) {
            *new_id = i;
            break;
        }
    }
    if (*new_id == -1) return NC_ENOMEM; /* out of memory */

    pnc_filelist[*new_id] = pncp;
    pnc_numfiles++;
    return NC_NOERR;
}

static void
del_from_PNCList(int ncid)
{
    /* validity of ncid should have been checked already */
    pnc_filelist[ncid] = NULL;
    pnc_numfiles--;
}

static PNC*
find_in_PNCList_by_name(const char* path)
{
    int i;
    PNC* pncp = NULL;
    for (i=0; i<NC_MAX_FILES; i++) {
         if (pnc_filelist[i] == NULL) continue;
         if (strcmp(pnc_filelist[i]->path, path) == 0) {
             pncp = pnc_filelist[i];
             break;
         }
    }
    return pncp;
}

int
PNC_check_id(int ncid, PNC **pncp)
{
    if (pnc_numfiles == 0 || ncid < 0 || ncid > NC_MAX_FILES)
        return NC_EBADID;

    *pncp = pnc_filelist[ncid];

    return NC_NOERR;
}

/*----< ncmpi_create() >-----------------------------------------------------*/
int
ncmpi_create(MPI_Comm    comm,
             const char *path,
             int         cmode,
             MPI_Info    info,
             int        *ncidp)
{
    int err;
    PNC *pncp;

    /* Use comde to tell the file format which is later used to select the
     * right dispatcher.
     */
    PNC_Dispatch *dispatcher = ncmpii_inq_dispatcher();

    /* check whether this path is already opened */
    pncp = find_in_PNCList_by_name(path);
    if (pncp != NULL) return NC_ENFILE;

    /* Create a PNC object and insert its dispatcher */
    pncp = (PNC*) malloc(sizeof(PNC));
    if (pncp == NULL) return NC_ENOMEM;
    pncp->path = (char*) malloc(strlen(path)+1);
    if (pncp->path == NULL) {
        free(pncp);
        return NC_ENOMEM;
    }
    strcpy(pncp->path, path);
    pncp->mode = cmode;
    pncp->dispatch = dispatcher;

    /* calling the create subroutine */
    err = dispatcher->create(comm, path, cmode, info, &pncp->ncp);
    if (err != NC_NOERR) {
        free(pncp->path);
        free(pncp);
        return err;
    }

    /* add to the PNCList and obtain ncid */
    err = add_to_PNCList(pncp, ncidp);
    if (err != NC_NOERR) {
        free(pncp->path);
        free(pncp);
        return err;
    }
    return NC_NOERR;
}

/*----< ncmpi_open() >-------------------------------------------------------*/
int
ncmpi_open(MPI_Comm    comm,
           const char *path,
           int         omode,
           MPI_Info    info,
           int        *ncidp)
{
    int err;
    PNC *pncp;

    /* Check the file signature to tell the file format which is later used to
     * select the right dispatcher.
     */
    PNC_Dispatch *dispatcher = ncmpii_inq_dispatcher();

    /* Create a PNC object and insert its dispatcher */
    pncp = (PNC*) malloc(sizeof(PNC));
    if (pncp == NULL) return NC_ENOMEM;
    pncp->path = (char*) malloc(strlen(path)+1);
    if (pncp->path == NULL) {
        free(pncp);
        return NC_ENOMEM;
    }
    strcpy(pncp->path, path);
    pncp->mode = omode;
    pncp->dispatch = dispatcher;

    /* calling the create subroutine */
    err = dispatcher->open(comm, path, omode, info, &pncp->ncp);
    if (err != NC_NOERR) {
        free(pncp->path);
        free(pncp);
        return err;
    }

    /* add to the PNCList and obtain ncid */
    err = add_to_PNCList(pncp, ncidp);
    if (err != NC_NOERR) {
        free(pncp->path);
        free(pncp);
        return err;
    }
    return NC_NOERR;
}

/*----< ncmpi_close() >------------------------------------------------------*/
int
ncmpi_close(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_close() */
    err = pncp->dispatch->close(pncp->ncp);
    if (err != NC_NOERR) return err;

    /* Remove from the PNCList */
    del_from_PNCList(ncid);

    /* free the PNC object */
    free(pncp->path);
    free(pncp);

    return NC_NOERR;
}

/*----< ncmpi_enddef() >-----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_enddef(int ncid) {
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_enddef() */
    err = pncp->dispatch->enddef(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi__enddef() >----------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi__enddef(int        ncid,
              MPI_Offset h_minfree,
              MPI_Offset v_align,
              MPI_Offset v_minfree,
              MPI_Offset r_align)
{   
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi__enddef() */
    err = pncp->dispatch->_enddef(pncp->ncp, h_minfree, v_align, v_minfree, r_align);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_redef() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_redef(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_redef() */
    err = pncp->dispatch->redef(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_sync() >-------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpi_sync(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync() */
    err = pncp->dispatch->sync(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_abort() >------------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_abort(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_abort() */
    err = pncp->dispatch->abort(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_set_fill() >---------------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_set_fill(int  ncid,
               int  fill_mode,
               int *old_fill_mode)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_set_fill() */
    err = pncp->dispatch->set_fill(pncp->ncp, fill_mode, old_fill_mode);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_format() >-------------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_format(int  ncid,
                 int *formatp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->format == 5) {
        *formatp = NC_FORMAT_CDF5;
    } else if (pncp->format == 2) {
        *formatp = NC_FORMAT_CDF2;
    } else if (pncp->format == 1) {
        *formatp = NC_FORMAT_CLASSIC;
    } else {
        /* this should not happen, because if ncid is valid, checking whether
         * the file is in a supported CDF format should have already been done
         * at ncmpi_open or ncmpi_create
         */
        *formatp = NC_FORMAT_UNKNOWN;
    }

    return NC_NOERR;
}

/*----< ncmpi_inq_file_format() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_format(const char *filename,
                      int        *formatp) /* out */
{
#ifdef _USE_NCMPI
    int ncid, err;
    NC *ncp;

    /* open file for reading its header */
    err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        if (err == NC_ENOTNC3)
            DEBUG_ASSIGN_ERROR(*formatp, NC_FORMAT_NETCDF4)
        else if (err == NC_ENOTNC)
            DEBUG_ASSIGN_ERROR(*formatp, NC_FORMAT_UNKNOWN)
        return err;
    }

    err = ncmpii_NC_check_id(ncid, &ncp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    if (ncp->format == 5) {
        *formatp = NC_FORMAT_CDF5;
    } else if (ncp->format == 2) {
        *formatp = NC_FORMAT_CDF2;
    } else {  /* if (ncp->format == 1) */
        *formatp = NC_FORMAT_CLASSIC;
    }
    err = ncmpi_close(ncid);

    return err;
#else
    char *cdf_signature="CDF";
    char *hdf5_signature="\211HDF\r\n\032\n";
    char signature[8];
    int fd;
    ssize_t rlen;

    *formatp = NC_FORMAT_UNKNOWN;

    if ((fd = open(filename, O_RDONLY, 0700)) == -1) {
             if (errno == ENOENT)       return NC_ENOENT;
        else if (errno == EACCES)       return NC_EACCESS;
        else if (errno == ENAMETOOLONG) return NC_EBAD_FILE;
        else                            return NC_EFILE;
    }
    /* get first 8 bytes of file */
    rlen = read(fd, signature, 8);
    if (rlen != 8) {
        close(fd); /* ignore error */
        return NC_EFILE;
    }
    if (close(fd) == -1) {
        return NC_EFILE;
    }

    if (memcmp(signature, hdf5_signature, 8) == 0) {
        /* whether the file is NC_FORMAT_NETCDF4_CLASSIC is determined by HDF5
         * attribute "_nc3_strict" which requires a call to H5Aget_name(). Here
         * we do not distinquish NC_CLASSIC_MODEL, but simply return NETCDF4
         * format.
         */
        *formatp = NC_FORMAT_NETCDF4;
    }
    else if (memcmp(signature, cdf_signature, 3) == 0) {
             if (signature[3] == 5)  *formatp = NC_FORMAT_CDF5;
        else if (signature[3] == 2)  *formatp = NC_FORMAT_CDF2;
        else if (signature[3] == 1)  *formatp = NC_FORMAT_CLASSIC;
    }

    return NC_NOERR;
#endif
}

/*----< ncmpi_inq() >------------------------------------------------------*/
int
ncmpi_inq(int  ncid,
          int *ndimsp,
          int *nvarsp,
          int *nattsp,
          int *xtendimp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_abort() */
    err = pncp->dispatch->inq(pncp->ncp, ndimsp, nvarsp, nattsp, xtendimp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_ndims() >--------------------------------------------------*/
int
ncmpi_inq_ndims(int  ncid,
                int *ndimsp)
{
    return ncmpi_inq(ncid, ndimsp, NULL, NULL, NULL);
}

/*----< ncmpi_inq_nvars() >--------------------------------------------------*/
int
ncmpi_inq_nvars(int  ncid,
                int *nvarsp)
{
    return ncmpi_inq(ncid, NULL, nvarsp, NULL, NULL);
}

/*----< ncmpi_inq_natts() >--------------------------------------------------*/
int
ncmpi_inq_natts(int  ncid,
                int *nattsp)
{
    return ncmpi_inq(ncid, NULL, NULL, nattsp, NULL);
}

/*----< ncmpi_inq_unlimdim() >-----------------------------------------------*/
int
ncmpi_inq_unlimdim(int  ncid,
                   int *unlimdimidp)
{
    return ncmpi_inq(ncid, NULL, NULL, NULL, unlimdimidp);
}

/*----< ncmpi_begin_indep_data() >-------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_begin_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_begin_indep_data() */
    err = pncp->dispatch->begin_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_end_indep_data() >---------------------------------------------*/
/* This is a collective subroutine. */
int
ncmpi_end_indep_data(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_end_indep_data() */
    err = pncp->dispatch->end_indep_data(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

