/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifndef _H_MACRO
#define _H_MACRO

#define WRITE_REQ 0
#define READ_REQ  1

#define INDEP_IO 0
#define COLL_IO  1

#ifndef MAX
#define MAX(mm,nn) (((mm) > (nn)) ? (mm) : (nn))
#endif
#ifndef MIN
#define MIN(mm,nn) (((mm) < (nn)) ? (mm) : (nn))
#endif

void *NCI_Malloc_fn(size_t size, int lineno, const char *fname);
void *NCI_Calloc_fn(size_t nelem, size_t elsize, int lineno, const char *fname);
void *NCI_Realloc_fn(void *ptr, size_t size, int lineno, const char *fname);
void NCI_Free_fn(void *ptr, int lineno, const char *fname);

#define NCI_Malloc(a) NCI_Malloc_fn(a,__LINE__,__FILE__)
#define NCI_Calloc(a,b) NCI_Calloc_fn(a,b,__LINE__,__FILE__)
#define NCI_Realloc(a,b) NCI_Realloc_fn(a,b,__LINE__,__FILE__)
#define NCI_Free(a) NCI_Free_fn(a,__LINE__,__FILE__)


#define CHECK_MPI_ERROR(str, err) {                                   \
    if (mpireturn != MPI_SUCCESS) {                                   \
        char errorString[MPI_MAX_ERROR_STRING];                       \
        int rank, errorStringLen;                                     \
        MPI_Comm_rank(ncp->nciop->comm, &rank);                       \
        MPI_Error_string(mpireturn, errorString, &errorStringLen);    \
        printf("%2d: MPI Failure at line %d of %s (%s : %s)\n",       \
               rank, __LINE__, __FILE__, str, errorString);           \
        return err;                                                   \
    }                                                                 \
}

/* API error will terminate the API call, not the entire program */
#define CHECK_NCID {                                          \
    status = ncmpii_NC_check_id(ncid, &ncp);                  \
    if (status != NC_NOERR) { /* API error */                 \
        /* uncomment to print debug message                   \
        int rank;                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                 \
        printf("%2d: Invalid ncid(%d) at line %d of %s\n",    \
               rank, ncid, __LINE__, __FILE__);               \
        */                                                    \
        return status;                                        \
    }                                                         \
}

#define CHECK_VARID(varid, varp) {                            \
    varp = ncmpii_NC_lookupvar(ncp, varid);                   \
    if (varp == NULL) { /* API error */                       \
        /* uncomment to print debug message                   \
        int rank;                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                 \
        printf("%2d: Invalid varid(%d) at line %d of %s\n",   \
               rank, varid, __LINE__, __FILE__);              \
        */                                                    \
        return NC_ENOTVAR;                                    \
    }                                                         \
}

#define GET_TOTAL_NUM_ELEMENTS {                    \
    int ndims = varp->ndims;                        \
    if (ndims == 0)                                 \
        nelems = 1;                                 \
    else if (!IS_RECVAR(varp))                      \
        nelems = varp->dsizes[0];                   \
    else if (ndims > 1)                             \
        nelems = ncp->numrecs * varp->dsizes[1];    \
    else                                            \
        nelems = ncp->numrecs;                      \
}

#define GET_NUM_ELEMENTS {           \
    int _i;                           \
    nelems = 1;                      \
    for (_i=0; _i<varp->ndims; _i++)    \
        nelems *= count[_i];          \
}

#define CHECK_WRITE_PERMISSION {                                   \
    if (NC_readonly(ncp)) { /* API error */                        \
        int rank;                                                  \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                      \
        printf("%2d: No file write permission at line %d of %s\n", \
               rank, __LINE__, __FILE__);                          \
        return NC_EPERM;                                           \
    }                                                              \
}

#define CHECK_DATATYPE(datatype, ptype, esize, cnelems, iscontig) {        \
    int isderived;                                                         \
    err = ncmpii_dtype_decode(datatype, &(ptype), &(esize), &(cnelems),    \
                              &isderived, &iscontig);                      \
    if (err != NC_NOERR) { /* API error */                                 \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
        printf("%2d: datatype decode error at line %d of %s\n",            \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        goto err_check;                                                    \
    }                                                                      \
}

#define CHECK_ECHAR(varp) {                                                \
    /* unable to type convert for char type */                             \
    if ( ncmpii_echar((varp)->type, ptype) ) { /* API error */             \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
        printf("%2d: datatype cannot convert to CHAR at line %d of %s\n",  \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        err = NC_ECHAR;                                                    \
        goto err_check;                                                    \
    }                                                                      \
}

#define CHECK_NELEMS(varp, cnelems, count, bufcount, nelems, nbytes) {     \
    int _i;                                                                \
    /* cnelems is calculated from the number of elements in datatype */    \
    cnelems *= bufcount;                                                   \
                                                                           \
    /* nelems is calculated from count[] */                                \
    nelems = 1;                                                            \
    for (_i=0; _i<(varp)->ndims; _i++) {                                   \
        if (count[_i] < 0) { /* API error */                               \
            err = NC_ENEGATIVECNT;                                         \
            goto err_check;                                                \
        }                                                                  \
        nelems *= count[_i];                                               \
    }                                                                      \
                                                                           \
    /* check mismatch between cnelems and nelems */                        \
    if (nelems != cnelems) {                                               \
        if (warning == NC_NOERR)                                           \
            warning = NC_EIOMISMATCH;                                      \
        (nelems>cnelems) ? (nelems=cnelems) : (cnelems=nelems);            \
    }                                                                      \
    /* now nelems == cnelems */                                            \
                                                                           \
    /* nbytes is the amount of this request in bytes */                    \
    nbytes = nelems * (varp)->xsz;                                         \
    if (nbytes < 0) { /* API error */                                      \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
        printf("%2d: Error - negative request amount at line %d of %s\n",  \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        err = NC_ENEGATIVECNT;                                             \
        goto err_check;                                                    \
    }                                                                      \
}

#define DATATYPE_GET_CONVERT(vartype, inbuf, outbuf, cnelems, ptype) {       \
    switch( vartype ) {                                                      \
        case NC_BYTE:                                                        \
            status = ncmpii_x_getn_schar(inbuf, outbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_SHORT:                                                       \
            status = ncmpii_x_getn_short(inbuf, outbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_INT:                                                         \
            status = ncmpii_x_getn_int(inbuf, outbuf, cnelems, ptype);       \
            break;                                                           \
        case NC_FLOAT:                                                       \
            status = ncmpii_x_getn_float(inbuf, outbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_DOUBLE:                                                      \
            status = ncmpii_x_getn_double(inbuf, outbuf, cnelems, ptype);    \
            break;                                                           \
        default:                                                             \
            break;                                                           \
    }                                                                        \
}

#define DATATYPE_PUT_CONVERT(vartype, outbuf, inbuf, cnelems, ptype) {       \
    switch( vartype ) {                                                      \
        case NC_BYTE:                                                        \
            status = ncmpii_x_putn_schar(outbuf, inbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_SHORT:                                                       \
            status = ncmpii_x_putn_short(outbuf, inbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_INT:                                                         \
            status = ncmpii_x_putn_int(outbuf, inbuf, cnelems, ptype);       \
            break;                                                           \
        case NC_FLOAT:                                                       \
            status = ncmpii_x_putn_float(outbuf, inbuf, cnelems, ptype);     \
            break;                                                           \
        case NC_DOUBLE:                                                      \
            status = ncmpii_x_putn_double(outbuf, inbuf, cnelems, ptype);    \
            break;                                                           \
        default:                                                             \
            break;                                                           \
    }                                                                        \
}

#define GET_FULL_DIMENSIONS {                                                \
    int _i;                                                                  \
    start = (MPI_Offset*) NCI_Malloc(2 * varp->ndims * sizeof(MPI_Offset));  \
    count = start + varp->ndims;                                             \
                                                                             \
    for (_i=0; _i<varp->ndims; _i++) {                                       \
        NC_dim *dimp;                                                        \
        dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t)varp->dimids[_i]); \
        if (dimp->size == NC_UNLIMITED)                                      \
            count[_i] = NC_get_numrecs(ncp);                                  \
        else                                                                 \
            count[_i] = dimp->size;                                           \
        start[_i] = 0;                                                        \
    }                                                                        \
}

#define GET_ONE_COUNT {                                                      \
    int _i;                                                                   \
    count = (MPI_Offset*) NCI_Malloc(varp->ndims * sizeof(MPI_Offset));      \
    for (_i=0; _i<varp->ndims; _i++)                                            \
        count[_i] = 1;                                                        \
}

#define CHECK_INDEP_FH {                                                      \
    /* check to see that the independent MPI file handle is opened */         \
    status =                                                                  \
    ncmpii_check_mpifh(ncp, &(ncp->nciop->independent_fh), MPI_COMM_SELF, 0); \
    if (status != NC_NOERR) {                                                 \
        /* uncomment to print debug message                                   \
        int rank;                                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                 \
        printf("%2d: Error - MPI indep file handle at line %d of %s\n",       \
               rank, __LINE__, __FILE__);                                     \
        */                                                                    \
        return status;                                                        \
    }                                                                         \
}

#define CHECK_COLLECTIVE_FH {                                                 \
    /* check to see that the collective MPI file handle is opened */          \
    status =                                                                  \
    ncmpii_check_mpifh(ncp, &(ncp->nciop->collective_fh),ncp->nciop->comm,1); \
    if (status != NC_NOERR) {                                                 \
        /* uncomment to print debug message                                   \
        int rank;                                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                 \
        printf("%2d: Error - MPI collective file handle at line %d of %s\n",  \
               rank, __LINE__, __FILE__);                                     \
        */                                                                    \
        return status;                                                        \
    }                                                                         \
} 

#endif
