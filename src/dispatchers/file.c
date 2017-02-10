
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
    int default_format, err;
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

    /* set the file format version based on the create mode, cmode */
    ncmpi_inq_default_format(&default_format);
    if (cmode & NC_64BIT_DATA) {
        pncp->format = 5;
    } else if (cmode & NC_64BIT_OFFSET) {
        pncp->format = 2;
    } else {
        if (default_format == NC_FORMAT_CDF5) {
            pncp->format = 5;
        }
        else if (default_format == NC_FORMAT_CDF2) {
            pncp->format = 2;
        }
        else {
            pncp->format = 1;
        }
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
    int format, err;
    PNC *pncp;

    /* Check the file signature to tell the file format which is later used to
     * select the right dispatcher.
     */
    err = ncmpi_inq_file_format(path, &format);
    if (err != NC_NOERR) return err;

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

    if (format == NC_FORMAT_CDF5) pncp->format = 5;
    else if (format == NC_FORMAT_CDF2) pncp->format = 2;
    else if (format == NC_FORMAT_CLASSIC) pncp->format = 1;
    else if (format == NC_FORMAT_NETCDF4) pncp->format = 4;

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
    } else if (pncp->format == 4) {
        *formatp = NC_FORMAT_NETCDF4;
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
    PNC *pncp;

    /* open file for reading its header */
    err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        if (err == NC_ENOTNC3)
            *formatp = NC_FORMAT_NETCDF4;
        else if (err == NC_ENOTNC)
            *formatp = NC_FORMAT_UNKNOWN;
        return err;
    }

    /* obtain pncp object pointer */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->format == 5) {
        *formatp = NC_FORMAT_CDF5;
    } else if (pncp->format == 2) {
        *formatp = NC_FORMAT_CDF2;
    } else {  /* if (pncp->format == 1) */
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

/*----< ncmpi_inq_striping() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_striping(int ncid, int *striping_size, int *striping_count)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_striping() */
    err = pncp->dispatch->inq_striping(pncp->ncp, striping_size, striping_count);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_num_rec_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_rec_vars(int ncid, int *num_rec_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_num_rec_vars() */
    err = pncp->dispatch->inq_num_rec_vars(pncp->ncp, num_rec_varsp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_num_fix_vars() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_num_fix_vars(int ncid, int *num_fix_varsp)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_num_fix_vars() */
    err = pncp->dispatch->inq_num_fix_vars(pncp->ncp, num_fix_varsp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_recsize() >-----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_recsize(int ncid, MPI_Offset *recsize)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_recsize() */
    err = pncp->dispatch->inq_recsize(pncp->ncp, recsize);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_version() >-----------------------------------------------*/
int
ncmpi_inq_version(int ncid, int *nc_mode)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    if (pncp->format == 5) {
        *nc_mode = NC_64BIT_DATA;
    } else if (pncp->format == 2) {
        *nc_mode = NC_64BIT_OFFSET;
    } else if (pncp->format == 1) {
        *nc_mode = NC_CLASSIC_MODEL;
    }

    return NC_NOERR;
}

/*----< ncmpi_inq_put_size() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_put_size(int ncid, MPI_Offset *size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_put_size() */
    err = pncp->dispatch->inq_put_size(pncp->ncp, size);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_get_size() >----------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_get_size(int ncid, MPI_Offset *size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_get_size() */
    err = pncp->dispatch->inq_get_size(pncp->ncp, size);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_header_size() >--------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_size(int ncid, MPI_Offset *size)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_header_size() */
    err = pncp->dispatch->inq_header_size(pncp->ncp, size);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_header_extent() >------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_header_extent(int ncid, MPI_Offset *extent)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_inq_header_extent() */
    err = pncp->dispatch->inq_header_extent(pncp->ncp, extent);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/*----< ncmpi_inq_file_info() >---------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_file_info(int ncid, MPI_Info *info)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_end_indep_data() */
    err = pncp->dispatch->inq_file_info(pncp->ncp, info);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

int
ncmpi_get_file_info(int ncid, MPI_Info *info)
{
    return ncmpi_inq_file_info(ncid, info);
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

/*----< ncmpi_sync_numrecs() >-----------------------------------------------*/
/* this API is collective, but can be called in independent data mode.
 * Note numrecs is always sync-ed in memory and update in file in collective
 * data mode.
 */
int
ncmpi_sync_numrecs(int ncid)
{
    int err;
    PNC *pncp;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_sync_numrecs() */
    err = pncp->dispatch->sync_numrecs(pncp->ncp);
    if (err != NC_NOERR) return err;

    return NC_NOERR;
}

/* The const string below is for the RCS ident(1) command to find a string like
 * "\044Id: \100(#) PnetCDF library version 1.4.0 of 16 Nov 2013 $"
 * in the library file (libpnetcdf.a).
 *
 * This string must be made a global variable. Otherwise, it won't work
 * when compiled with optimization options, e.g. -O2
 */
char const pnetcdf_libvers[] =
        "\044Id: \100(#) PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE" $";

/* a cleaner version for running command "strings", e.g.
 * % strings libpnetcdf.a | grep "PnetCDF library version"
 * or
 * % strings a.out | grep "PnetCDF library version"
 */
char pnetcdf_lib_vers[] = "PnetCDF library version "PNETCDF_VERSION" of "PNETCDF_RELEASE_DATE;

/* pnetcdf_libvers is slightly different from the one returned from
 * ncmpi_inq_libvers(). The string pnetcdf_libvers is for command "ident" to
 * use. People can run command ident libpnetcdf.a to obtain the version of a
 * library (or an executable built from that library). In PnetCDF case, the
 * command will print the string of pnetcdf_libvers. Command "ident' looks for
 * a specific keyword pattern and print it. See man page of ident.
 *
 * The API ncmpi_inq_libvers() below on the other hand returns a string to be
 * used by the utility tools like ncmpidump, ncmpigen, etc. Check the last line
 * of output from command "ncmpidump -v".
 */

/*----< ncmpi_inq_libvers() >------------------------------------------------*/
const char*
ncmpi_inq_libvers(void) {

    /* match the style used by netCDF API nc_inq_libvers()
     * for example, "4.3.0 of Jun 16 2013 12:11:30 $" */
    /* we need some silly operation so the compiler will emit the otherwise
     * unused pnetcdf_libvers */
    if ((void *)pnetcdf_libvers != (void *)ncmpi_inq_libvers) {
	; /* do nothing */
    }
    return PNETCDF_VERSION " of " PNETCDF_RELEASE_DATE;
}

/* This is the default create format for ncmpi_create and nc__create.
 * The use of this file scope variable is not thread-safe.
 */
static int ncmpi_default_create_format = NC_FORMAT_CLASSIC;

/*----< ncmpi_set_default_format() >-----------------------------------------*/
/* This function sets a default create file format.
 * Valid formats are NC_FORMAT_CLASSIC, NC_FORMAT_CDF2, and NC_FORMAT_CDF5
 * This API is NOT collective, as there is no way to check against an MPI
 * communicator. It should be called by all MPI processes that intend to
 * create a file later. Consistency check will have to be done in other APIs.
 */
int
ncmpi_set_default_format(int format, int *old_formatp)
{
    /* Return existing format if desired. */
    if (old_formatp != NULL)
        *old_formatp = ncmpi_default_create_format;

    /* Make sure only valid format is set. */
    if (format != NC_FORMAT_CLASSIC &&
        format != NC_FORMAT_CDF2 &&
        format != NC_FORMAT_CDF5) {
        return NC_EINVAL;
    }
    ncmpi_default_create_format = format;

    return NC_NOERR;
}

/* returns a value suitable for a create flag.  Will return one or more of the
 * following values OR-ed together:
 * NC_64BIT_OFFSET, NC_CLOBBER, NC_LOCK, NC_SHARE */
int
ncmpi_inq_default_format(int *formatp)
{
    if (formatp == NULL) return NC_EINVAL;

    *formatp = ncmpi_default_create_format;
    return NC_NOERR;
}

/*----< ncmpi_inq_files_opened() >-------------------------------------------*/
/* This is an independent subroutine. */
int
ncmpi_inq_files_opened(int *num, int *ncids)
{
    int i;

    if (num == NULL) return NC_EINVAL;

    *num = 0;
    for (i=0; i<NC_MAX_FILES; i++) {
        if (pnc_filelist[i] != NULL) {
            if (ncids != NULL) /* ncids can be NULL */
                ncids[*num] = i;
            (*num)++;
        }
    }
    return NC_NOERR;
}

