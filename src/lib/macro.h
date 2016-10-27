/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _H_MACRO
#define _H_MACRO

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

void *NCI_Malloc_fn(size_t size, int lineno, const char *func, const char *fname);
void *NCI_Calloc_fn(size_t nelem, size_t elsize, int lineno, const char *func, const char *fname);
void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *func, const char *fname);
void  NCI_Free_fn(void *ptr, int lineno, const char *func, const char *fname);

#define NCI_Malloc(a)    NCI_Malloc_fn(a,__LINE__,__func__,__FILE__)
#define NCI_Calloc(a,b)  NCI_Calloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Realloc(a,b) NCI_Realloc_fn(a,b,__LINE__,__func__,__FILE__)
#define NCI_Free(a)      NCI_Free_fn(a,__LINE__,__func__,__FILE__)


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

#ifdef PNC_DEBUG
#define DEBUG_RETURN_ERROR(err) {                                       \
    char *env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");               \
    if (env_str != NULL && *env_str != '0') {                           \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpii_err_code_name(err),__LINE__,__func__,__FILE__);    \
    }                                                                   \
    return err;                                                         \
}
#define DEBUG_ASSIGN_ERROR(status, err) {                               \
    char *env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");               \
    if (env_str != NULL && *env_str != '0') {                           \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpii_err_code_name(err),__LINE__,__func__,__FILE__);    \
    }                                                                   \
    status = err;                                                       \
}
#define DEBUG_TRACE_ERROR {                                             \
    char *env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");               \
    if (env_str != NULL && *env_str != '0') {                           \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpii_err_code_name(err),__LINE__,__func__,__FILE__);    \
    }                                                                   \
}
#else
#define DEBUG_RETURN_ERROR(err) return err;
#define DEBUG_ASSIGN_ERROR(status, err) status = err;
#define DEBUG_TRACE_ERROR
#endif

#define GET_ONE_COUNT(count) {                                                \
    int _i;                                                                   \
    count = (MPI_Offset*) NCI_Malloc((size_t)varp->ndims * SIZEOF_MPI_OFFSET);\
    for (_i=0; _i<varp->ndims; _i++)                                          \
        count[_i] = 1;                                                        \
}

#ifdef ENABLE_SUBFILING
#define GET_FULL_DIMENSIONS(start, count) {                                   \
    int _i;                                                                   \
    int _ndims = (varp->num_subfiles>1?varp->ndims_org:varp->ndims);          \
    start = (MPI_Offset*) NCI_Malloc((size_t)_ndims*2*SIZEOF_MPI_OFFSET);     \
    count = start + _ndims;                                                   \
                                                                              \
    for (_i=0; _i<_ndims; _i++) {                                             \
        NC_dim *dimp;                                                         \
        int _dimid;                                                           \
        _dimid = (varp->num_subfiles>1)?varp->dimids_org[_i]:varp->dimids[_i];\
        dimp = ncmpii_elem_NC_dimarray(&ncp->dims, _dimid);                   \
        if (dimp->size == NC_UNLIMITED)                                       \
            count[_i] = NC_get_numrecs(ncp);                                  \
        else                                                                  \
            count[_i] = dimp->size;                                           \
        start[_i] = 0;                                                        \
    }                                                                         \
}
#else /* without subfiling */
#define GET_FULL_DIMENSIONS(start, count) {                                   \
    int _i=0;                                                                 \
    start = (MPI_Offset*) NCI_Malloc((size_t)varp->ndims*2*SIZEOF_MPI_OFFSET);\
    count = start + varp->ndims;                                              \
                                                                              \
    if (IS_RECVAR(varp)) { /* find current numrec if varp is record var */    \
        count[0] = NC_get_numrecs(ncp);                                       \
        start[0] = 0;                                                         \
        _i = 1;                                                               \
    }                                                                         \
    for (; _i<varp->ndims; _i++) {                                            \
        count[_i] = varp->shape[_i];                                          \
        start[_i] = 0;                                                        \
    }                                                                         \
}
#endif

#define DATATYPE_GET_CONVERT(cdf_ver,xtype,xbuf,ibuf,bnelems,itype,err) {     \
    /* xtype is the NC variable's external data type stored in the nc file    \
     * itype (internal type) is the user I/O buffer data type (MPI_Datatype)  \
     * xbuf is the buffer containing the data read from the nc file in the    \
     * external representation, to be converted to internal representation r  \
     * into user buffer, ibuf */                                              \
    switch(xtype) {                                                           \
        case NC_BYTE:                                                         \
            err = ncmpii_x_getn_NC_BYTE(cdf_ver,xbuf,ibuf,bnelems,itype);     \
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_x_getn_NC_UBYTE(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_x_getn_NC_SHORT(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_x_getn_NC_USHORT(xbuf,ibuf,bnelems,itype);           \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_x_getn_NC_INT(xbuf,ibuf,bnelems,itype);              \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_x_getn_NC_UINT(xbuf,ibuf,bnelems,itype);             \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_x_getn_NC_FLOAT(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_x_getn_NC_DOUBLE(xbuf,ibuf,bnelems,itype);           \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_x_getn_NC_INT64(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_x_getn_NC_UINT64(xbuf,ibuf,bnelems,itype);           \
            break;                                                            \
        default:                                                              \
            err = NC_EBADTYPE;                                                \
            break;                                                            \
    }                                                                         \
}

#define DATATYPE_PUT_CONVERT(cdf_ver,xtype,xbuf,ibuf,cnelems,itype,fillp,err) {\
    /* xtype is the NC variable's external data type stored in the nc file    \
     * itype (internal type)is the user I/O buffer data type (MPI_Datatype)   \
     * ibuf is user buffer containing the data in internal format of itype,   \
     * they meed to be converted to external format to external buffer, xbuf, \
     * so xbuf can be used to write to the nc file through the call to        \
     * MPI_FIle_write call */                                                 \
    switch(xtype) {                                                           \
        case NC_BYTE:                                                         \
            err = ncmpii_x_putn_NC_BYTE(cdf_ver,xbuf,ibuf,cnelems,itype,fillp);\
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_x_putn_NC_UBYTE(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_x_putn_NC_SHORT(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_x_putn_NC_USHORT(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_x_putn_NC_INT(xbuf,ibuf,cnelems,itype,fillp);        \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_x_putn_NC_UINT(xbuf,ibuf,cnelems,itype,fillp);       \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_x_putn_NC_FLOAT(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_x_putn_NC_DOUBLE(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_x_putn_NC_INT64(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_x_putn_NC_UINT64(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        default:                                                              \
            err = NC_EBADTYPE;                                                \
            break;                                                            \
    }                                                                         \
}

#endif
