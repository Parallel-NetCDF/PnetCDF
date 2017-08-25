/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _MACRO_H
#define _MACRO_H

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
        dimp = ncp->dims.value[_dimid];                                       \
        if (dimp->size == NC_UNLIMITED)                                       \
            count[_i] = ncp->numrecs;                                         \
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
        count[0] = ncp->numrecs;                                              \
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
            err = ncmpii_getn_NC_BYTE(cdf_ver,xbuf,ibuf,bnelems,itype);     \
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_getn_NC_UBYTE(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_getn_NC_SHORT(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_getn_NC_USHORT(xbuf,ibuf,bnelems,itype);           \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_getn_NC_INT(xbuf,ibuf,bnelems,itype);              \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_getn_NC_UINT(xbuf,ibuf,bnelems,itype);             \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_getn_NC_FLOAT(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_getn_NC_DOUBLE(xbuf,ibuf,bnelems,itype);           \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_getn_NC_INT64(xbuf,ibuf,bnelems,itype);            \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_getn_NC_UINT64(xbuf,ibuf,bnelems,itype);           \
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
            err = ncmpii_putn_NC_BYTE(cdf_ver,xbuf,ibuf,cnelems,itype,fillp);\
            break;                                                            \
        case NC_UBYTE:                                                        \
            err = ncmpii_putn_NC_UBYTE(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_SHORT:                                                        \
            err = ncmpii_putn_NC_SHORT(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_USHORT:                                                       \
            err = ncmpii_putn_NC_USHORT(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        case NC_INT:                                                          \
            err = ncmpii_putn_NC_INT(xbuf,ibuf,cnelems,itype,fillp);        \
            break;                                                            \
        case NC_UINT:                                                         \
            err = ncmpii_putn_NC_UINT(xbuf,ibuf,cnelems,itype,fillp);       \
            break;                                                            \
        case NC_FLOAT:                                                        \
            err = ncmpii_putn_NC_FLOAT(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_DOUBLE:                                                       \
            err = ncmpii_putn_NC_DOUBLE(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        case NC_INT64:                                                        \
            err = ncmpii_putn_NC_INT64(xbuf,ibuf,cnelems,itype,fillp);      \
            break;                                                            \
        case NC_UINT64:                                                       \
            err = ncmpii_putn_NC_UINT64(xbuf,ibuf,cnelems,itype,fillp);     \
            break;                                                            \
        default:                                                              \
            err = NC_EBADTYPE;                                                \
            break;                                                            \
    }                                                                         \
}

#endif
