/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _H_MACRO
#define _H_MACRO

#define WRITE_REQ 0
#define READ_REQ  1

#define INDEP_IO 0
#define COLL_IO  1
#define INDEP_COLL_IO  -1

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

void *NCI_Malloc_fn(size_t size, int lineno, const char *fname);
void *NCI_Calloc_fn(size_t nelem, size_t elsize, int lineno, const char *fname);
void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *fname);
void  NCI_Free_fn(void *ptr, int lineno, const char *fname);

#define NCI_Malloc(a)    NCI_Malloc_fn(a,__LINE__,__FILE__)
#define NCI_Calloc(a,b)  NCI_Calloc_fn(a,b,__LINE__,__FILE__)
#define NCI_Realloc(a,b) NCI_Realloc_fn(a,b,__LINE__,__FILE__)
#define NCI_Free(a)      NCI_Free_fn(a,__LINE__,__FILE__)


#define CHECK_MPI_ERROR(mpi_errorcode, err_msg, nc_err) {                     \
    if (mpi_errorcode != MPI_SUCCESS) {                                       \
        char errorString[MPI_MAX_ERROR_STRING];                               \
        int rank, errorStringLen;                                             \
        MPI_Comm_rank(ncp->nciop->comm, &rank);                               \
        MPI_Error_string(mpi_errorcode, errorString, &errorStringLen);        \
        printf("%2d: MPI Failure at line %d of %s (%s : %s)\n",               \
               rank, __LINE__, __FILE__, err_msg, errorString);               \
        mpi_err = nc_err;                                                     \
    }                                                                         \
}


#define SANITY_CHECK(ncid, ncp, varp, rw_flag, io_method, status) {            \
    int min_st;                                                                \
                                                                               \
    /* check if ncid is valid */                                               \
    status = ncmpii_NC_check_id(ncid, &ncp);                                   \
    if (status != NC_NOERR)                                                    \
        /* must return error now, collective APIs might hang if error occurs   \
           only on a subset of processes */                                    \
        return status;                                                         \
                                                                               \
    /* check if it is in define mode */                                        \
    if (NC_indef(ncp)) status = NC_EINDEFINE;                                  \
                                                                               \
    /* check if varid is valid */                                              \
    if (status == NC_NOERR) {                                                  \
        varp = ncmpii_NC_lookupvar(ncp, varid);                                \
        if (varp == NULL) status = NC_ENOTVAR;                                 \
    }                                                                          \
    /* check file write permission if this is write request */                 \
    if (status == NC_NOERR) {                                                  \
        if (rw_flag == WRITE_REQ && NC_readonly(ncp)) status = NC_EPERM;       \
    }                                                                          \
    /* check whether collective or independent mode */                         \
    if (status == NC_NOERR) {                                                  \
        if (io_method == INDEP_IO)                                             \
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->independent_fh),    \
                                        MPI_COMM_SELF, 0);                     \
        else if (io_method == COLL_IO)                                         \
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->collective_fh),     \
                                        ncp->nciop->comm, 1);                  \
        /* else if (io_method == INDEP_COLL_IO) */                             \
    }                                                                          \
    if (ncp->safe_mode == 1 && io_method == COLL_IO)                           \
        MPI_Allreduce(&status, &min_st, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);\
    else                                                                       \
        min_st = status;                                                       \
                                                                               \
    if (min_st != NC_NOERR)                                                    \
        return status;                                                         \
}

#define GET_ONE_COUNT(count) {                                                \
    int _i;                                                                   \
    count = (MPI_Offset*) NCI_Malloc(varp->ndims * sizeof(MPI_Offset));       \
    for (_i=0; _i<varp->ndims; _i++)                                          \
        count[_i] = 1;                                                        \
}

#define GET_TOTAL_NUM_ELEMENTS(nelems) {                                      \
    int ndims = varp->ndims;                                                  \
    if (ndims == 0)                                                           \
        nelems = 1;                                                           \
    else if (!IS_RECVAR(varp))                                                \
        nelems = varp->dsizes[0];                                             \
    else if (ndims > 1)                                                       \
        nelems = ncp->numrecs * varp->dsizes[1];                              \
    else                                                                      \
        nelems = ncp->numrecs;                                                \
}

#define GET_NUM_ELEMENTS(nelems) {                                            \
    int _i;                                                                   \
    nelems = 1;                                                               \
    for (_i=0; _i<varp->ndims; _i++)                                          \
        nelems *= count[_i];                                                  \
}

#define GET_FULL_DIMENSIONS(start, count) {                                   \
    int _i;                                                                   \
    start = (MPI_Offset*) NCI_Malloc(2 * varp->ndims * sizeof(MPI_Offset));   \
    count = start + varp->ndims;                                              \
                                                                              \
    for (_i=0; _i<varp->ndims; _i++) {                                        \
        NC_dim *dimp;                                                         \
        dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t)varp->dimids[_i]); \
        if (dimp->size == NC_UNLIMITED)                                       \
            count[_i] = NC_get_numrecs(ncp);                                  \
        else                                                                  \
            count[_i] = dimp->size;                                           \
        start[_i] = 0;                                                        \
    }                                                                         \
}

#define CHECK_NELEMS(varp, fnelems, fcount, bnelems, bufcount, nbytes, err) { \
    int _i;                                                                   \
    /* bnelems is the number of elements in the I/O buffer and the element    \
     * is of MPI primitive type. When input, bnelems is the number of         \
     * MPI primitive type elements in the user's derived data type.           \
     * Here, we make it the total MPI primitive type elements in the          \
     * user's buffer                                                          \
     */                                                                       \
    bnelems *= bufcount;                                                      \
                                                                              \
    /* fnelems is the total number of nc_type elements calculated from        \
     * fcount[]. fcount[] is the access count to the variable defined in      \
     * the netCDF file.                                                       \
     */                                                                       \
    fnelems = 1;                                                              \
    for (_i=0; _i<(varp)->ndims; _i++) {                                      \
        if (fcount[_i] < 0) { /* API error */                                 \
            err = NC_ENEGATIVECNT;                                            \
            goto err_check;                                                   \
            /* for collective call must return collectively */                \
        }                                                                     \
        fnelems *= fcount[_i];                                                \
    }                                                                         \
                                                                              \
    /* check mismatch between bnelems and fnelems */                          \
    if (fnelems != bnelems) {                                                 \
        if (warning == NC_NOERR)                                              \
            warning = NC_EIOMISMATCH;                                         \
        (fnelems>bnelems) ? (fnelems=bnelems) : (bnelems=fnelems);            \
        /* only handle partial of the request, smaller number of the two */   \
    }                                                                         \
    /* now fnelems == bnelems */                                              \
                                                                              \
    /* nbytes is the amount of this request in bytes */                       \
    nbytes = fnelems * (varp)->xsz;                                           \
    if (nbytes < 0) { /* API error */                                         \
        /* uncomment to print debug message                                   \
        int rank;                                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                 \
        printf("%2d: Error - negative request amount at line %d of %s\n",     \
               rank, __LINE__, __FILE__);                                     \
        */                                                                    \
        err = NC_ENEGATIVECNT;                                                \
        goto err_check;                                                       \
        /* cannot return now, for collective call must return collectively */ \
    }                                                                         \
}

#define DATATYPE_GET_CONVERT(vartype, inbuf, outbuf, cnelems, memtype, err) { \
    /* vartype is the variable's data type defined in the nc file             \
     * memtype is the I/O buffers data type (MPI_Datatype)  */                \
    switch(vartype) {                                                         \
        case NC_BYTE:                                                         \
            err = ncmpii_x_getn_schar(inbuf, outbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_x_getn_uchar(inbuf, outbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_x_getn_short(inbuf, outbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_x_getn_ushort(inbuf, outbuf, cnelems, memtype);      \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_x_getn_int(inbuf, outbuf, cnelems, memtype);         \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_x_getn_uint(inbuf, outbuf, cnelems, memtype);        \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_x_getn_float(inbuf, outbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_x_getn_double(inbuf, outbuf, cnelems, memtype);      \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_x_getn_int64(inbuf, outbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_x_getn_uint64(inbuf, outbuf, cnelems, memtype);      \
            break;                                                            \
        default:                                                              \
            err = NC_EBADTYPE;                                                \
            break;                                                            \
    }                                                                         \
}

#define DATATYPE_PUT_CONVERT(vartype, outbuf, inbuf, cnelems, memtype, err) { \
    /* vartype is the variable's data type defined in the nc file             \
     * memtype is the I/O buffers data type (MPI_Datatype)  */                \
    switch(vartype) {                                                         \
        case NC_BYTE:                                                         \
            err = ncmpii_x_putn_schar(outbuf, inbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_x_putn_uchar(outbuf, inbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_x_putn_short(outbuf, inbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_x_putn_ushort(outbuf, inbuf, cnelems, memtype);      \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_x_putn_int(outbuf, inbuf, cnelems, memtype);         \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_x_putn_uint(outbuf, inbuf, cnelems, memtype);        \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_x_putn_float(outbuf, inbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_x_putn_double(outbuf, inbuf, cnelems, memtype);      \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_x_putn_int64(outbuf, inbuf, cnelems, memtype);       \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_x_putn_uint64(outbuf, inbuf, cnelems, memtype);      \
            break;                                                            \
        default:                                                              \
            err = NC_EBADTYPE;                                                \
            break;                                                            \
    }                                                                         \
}

#endif
