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

#define CHECK_NCID {                                          \
    status = ncmpii_NC_check_id(ncid, &ncp);                  \
    if (status != NC_NOERR) {                                 \
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
    if (varp == NULL) {                                       \
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
    int i;                           \
    nelems = 1;                      \
    for (i=0; i<varp->ndims; i++)    \
        nelems *= count[i];          \
}

#define CHECK_WRITE_PERMISSION {                                   \
    if (NC_readonly(ncp)) {                                        \
        int rank;                                                  \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                      \
	printf("%2d: No file write permission at line %d of %s\n", \
               rank, __LINE__, __FILE__);                          \
        return NC_EPERM;                                           \
    }                                                              \
}

#define CHECK_DATATYPE(datatype, ptype, esize, cnelems, iscontig) {        \
    int isderived;                                                         \
    status = ncmpii_dtype_decode(datatype, &(ptype), &(esize), &(cnelems), \
                                 &isderived, &iscontig);                   \
    if (status != NC_NOERR) {                                              \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
	printf("%2d: datatype decode error at line %d of %s\n",            \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        return status;                                                     \
    }                                                                      \
}

#define CHECK_ECHAR(varp) {                                                \
    /* unable to type convert for char type */                             \
    if ( ncmpii_echar((varp)->type, ptype) ) {                             \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
	printf("%2d: datatype cannot convert to CHAR at line %d of %s\n",  \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        return NC_ECHAR;                                                   \
    }                                                                      \
}

#define CHECK_NELEMS(varp, cnelems, count, bufcount, nelems, nbytes) {     \
    /* cnelems is calculated from the number of elements in datatype */    \
    cnelems *= bufcount;                                                   \
                                                                           \
    /* nelems is calculated from count[] */                                \
    nelems = 1;                                                            \
    for (i=0; i<(varp)->ndims; i++) {                                      \
        if (count[i] < 0)                                                  \
            return NC_ENEGATIVECNT;                                        \
        nelems *= count[i];                                                \
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
    if (nbytes < 0) {                                                      \
        /* uncomment to print debug message                                \
        int rank;                                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                              \
	printf("%2d: Error - negative request amount at line %d of %s\n",  \
               rank, __LINE__, __FILE__);                                  \
        */                                                                 \
        return NC_ENEGATIVECNT;                                            \
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
    int i;                                                                   \
    start = (MPI_Offset*) malloc(2 * varp->ndims * sizeof(MPI_Offset));      \
    count = start + varp->ndims;                                             \
                                                                             \
    for (i=0; i<varp->ndims; i++) {                                          \
        NC_dim *dimp;                                                        \
        dimp = ncmpii_elem_NC_dimarray(&ncp->dims, (size_t)varp->dimids[i]); \
        if (dimp->size == NC_UNLIMITED)                                      \
            count[i] = NC_get_numrecs(ncp);                                  \
        else                                                                 \
            count[i] = dimp->size;                                           \
        start[i] = 0;                                                        \
    }                                                                        \
}

#define GET_ONE_COUNT {                                             \
    int i;                                                          \
    count = (MPI_Offset*) malloc(varp->ndims * sizeof(MPI_Offset)); \
    for (i=0; i<varp->ndims; i++)                                   \
        count[i] = 1;                                               \
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
