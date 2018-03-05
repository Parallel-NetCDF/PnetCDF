/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* Starting from version 1.8.0, PnetCDF checks file metadata consistency only
 * when safe mode is enabled and the checking is no longer done at enddef().
 * The check has been moved to where the new metadata is added, e.g.
 * ncmpi_def_dim, ncmpi_def_var, ncmpi_put_att, etc. Thus, below legacy codes
 * are for checking the entire metadata all at once.
 */
#ifdef _CHECK_HEADER_CONSISTENCY

#ifdef SIZEOF_INT
# if SIZEOF_INT == 4
#  define lld(x) (x)
# elif  SIZEOF_INT == 8
#  define lld(x) (long long)(x)
# endif
#endif

#define WARN_STR "Warning (inconsistent metadata):"

/*----< compare_dims() >-----------------------------------------------------*/
/* compare the local copy of dim_list against root's
 * If inconsistency is detected, overwrite local's with root's
 * this function is collective.
 */
static int
compare_dims(int          safe_mode,
             NC_dimarray *root_dim,
             NC_dimarray *local_dim)
{
    int i, err, status=NC_NOERR;

    if (root_dim->ndefined != local_dim->ndefined) {
        if (safe_mode)
            printf("%s number of dimensions (local=%d, root=%d)\n",
                   WARN_STR, local_dim->ndefined, root_dim->ndefined);
        DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_DIM_NUM)
    }

    for (i=0; i<root_dim->ndefined; i++) {

        if (i >= local_dim->ndefined) { /* if local list is shorter */
            /* copy root's dim to local */
            NC_dim *new_dim = ncmpio_dup_NC_dim(root_dim->value[i]);
            err = ncmpio_incr_NC_dimarray(local_dim, new_dim);
            if (err != NC_NOERR) return err; /* this is a fatal error */
            continue;
        }

        /* check dimension name */
        NC_string *root_name, *local_name;
         root_name =  root_dim->value[i]->name;
        local_name = local_dim->value[i]->name;

        err = NC_NOERR;
        if (root_name->nchars != local_name->nchars) {
            if (safe_mode)
                printf("%s dimension name length (local=%lld, root=%lld)\n",
                       WARN_STR, local_name->nchars, root_name->nchars);
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        }
        else if (memcmp(root_name->cp, local_name->cp, (size_t)root_name->nchars) != 0) {
            if (safe_mode)
                printf("%s dimension name (local=%s, root=%s)\n",
                       WARN_STR, local_name->cp, root_name->cp);
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_NAME)
        }
        else if (root_dim->value[i]->size != local_dim->value[i]->size) {
            /* check dimension size */
            if (safe_mode)
                printf("%s dimension %s's size (local=%lld, root=%lld)\n",
                       WARN_STR, root_dim->value[i]->name->cp,
                       root_dim->value[i]->size, local_dim->value[i]->size);
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_DIM_SIZE)
        }
        if (status == NC_NOERR) status = err;

        /* overwrite local's dim with root's */
        if (err != NC_NOERR) {
            ncmpio_free_NC_dim(local_dim->value[i]);
            local_dim->value[i] = ncmpio_dup_NC_dim(root_dim->value[i]);
        }
    }

    /* delete extra dimensions defined only in local copy */
    for (; i<local_dim->ndefined; i++)
        ncmpio_free_NC_dim(local_dim->value[i]);

    local_dim->ndefined = root_dim->ndefined;

#ifndef SEARCH_NAME_LINEARLY
    if (status != NC_NOERR) {
        /* dims are not consistent, must rebuild dim name lookup table */
        for (i=0; i<HASH_TABLE_SIZE; i++) {
            if (local_dim->nameT[i].num > 0)
                NCI_Free(local_dim->nameT[i].list);
            local_dim->nameT[i].num = 0;
            local_dim->nameT[i].list = NULL;
        }

        /* populate dim name lookup table */
        for (i=0; i<local_dim->ndefined; i++) {
            /* hash the dim name into a key for name lookup */
            int key = HASH_FUNC(local_dim->value[i]->name->cp);
            NC_nametable *nameT = &local_dim->nameT[key];
            if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
                nameT->list = (int*) NCI_Realloc(nameT->list,
                              (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
            nameT->list[nameT->num] = i;
            nameT->num++;
        }
    }
#endif

    return status;
}

/*----< compare_attrs() >-----------------------------------------------------*/
/* compare the local copy of attr_list against root's
 * If inconsistency is detected, overwrite local's with root's
 */
static int
compare_attrs(int           safe_mode,
              NC_attrarray *root_attr,
              NC_attrarray *local_attr)
{
    int i, j, err, status=NC_NOERR;
    char *msg;

    /* check if the numbers of attributes are the same */
    if (root_attr->ndefined != local_attr->ndefined) {
        if (safe_mode)
            printf("%s number of attributes (root=%d, local=%d)\n",
                   WARN_STR, root_attr->ndefined, local_attr->ndefined);
        DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_ATTR_NUM)
    }

    for (i=0; i<root_attr->ndefined; i++) {

        if (i >= local_attr->ndefined) { /* if local list is shorter */
            /* copy root's attr to local */
            NC_attr *new_attr = ncmpio_dup_NC_attr(root_attr->value[i]);
            err = ncmpio_incr_NC_attrarray(local_attr, new_attr);
            if (err != NC_NOERR) return err; /* a fatal error */
            continue;
        }

        NC_attr *v1 = root_attr->value[i];
        NC_attr *v2 = local_attr->value[i];
        char *name = v1->name->cp;

#define ATTR_WARN(msg, attr, root, local) \
    if (safe_mode) printf(msg, WARN_STR, attr, root, local);

#define ATTR_WARN_J(msg, attr, j, root, local) \
    if (safe_mode) printf(msg, WARN_STR, attr, j, root, local);

        err = NC_NOERR;
        if (v1->name->nchars != v2->name->nchars ||
            memcmp(name, v2->name->cp, (size_t)v1->name->nchars) != 0) {
            msg ="%s attribute %s (root=%s, local=%s)\n";
            ATTR_WARN(msg, "name", name, v2->name->cp)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
        }
        else if (v1->type != v2->type) {
            msg = "%s attribute \"%s\" type (root=%d, local=%d)\n";
            ATTR_WARN(msg, name, v1->type, v2->type)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_TYPE)
        }
        else if (v1->nelems != v2->nelems) {
            msg = "%s attribute \"%s\" length (root=%lld, local=%lld)\n";
            ATTR_WARN(msg, name, lld(v1->nelems), lld(v2->nelems))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_LEN)
        }
        else if (v1->xsz != v2->xsz) { /* internal check */
            msg = "%s attribute \"%s\" size (root=%lld, local=%lld)\n";
            ATTR_WARN(msg, name, lld(v1->xsz), lld(v2->xsz))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_SIZE)
        }
        /* hereinafter, we have v1->nelems == v2->nelems */
        else if (v1->type == NC_CHAR) {
            if (memcmp(v1->xvalue, v2->xvalue, (size_t)v1->nelems)) {
                msg = "%s attribute \"%s\" CHAR (root=%s, local=%s)\n";
                ATTR_WARN(msg, name, (char*)v1->xvalue, (char*)v2->xvalue);
                DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
            }
        }
        else if (v1->type == NC_BYTE) {
            schar *sba = (schar*) v1->xvalue;
            schar *sbb = (schar*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (sba[j] != sbb[j]) {
                    msg = "%s attribute \"%s\"[%d] BYTE (root=%hhdb, local=%hhdb)\n";
                    ATTR_WARN_J(msg, name, j, sba[j], sbb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_UBYTE) {
            uchar *uba = (uchar*) v1->xvalue;
            uchar *ubb = (uchar*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (uba[j] != ubb[j]) {
                    msg = "%s attribute \"%s\"[%d] UBYTE (root=%hhuub, local=%hhuub)\n";
                    ATTR_WARN_J(msg, name, j, uba[j], ubb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_SHORT) {
            short *ssa = (short*) v1->xvalue;
            short *ssb = (short*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (ssa[j] != ssb[j]) {
                    msg = "%s attribute \"%s\"[%d] SHORT (root=%hds, local=%hds)\n";
                    ATTR_WARN_J(msg, name, j, ssa[j], ssb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_USHORT) {
            ushort *usa = (ushort*) v1->xvalue;
            ushort *usb = (ushort*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (usa[j] != usb[j]) {
                    msg = "%s attribute \"%s\"[%d] USHORT (root=%huus, local=%huus)\n";
                    ATTR_WARN_J(msg, name, j, usa[j], usb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_INT) {
            int *sia = (int*) v1->xvalue;
            int *sib = (int*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (sia[j] != sib[j]) {
                    msg = "%s attribute \"%s\"[%d] INT (root=%d, local=%d)\n";
                    ATTR_WARN_J(msg, name, j, sia[j], sib[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_UINT) {
            uint *uia = (uint*) v1->xvalue;
            uint *uib = (uint*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (uia[j] != uib[j]) {
                    msg = "%s attribute \"%s\"[%d] UINT (root=%uu, local=%uu)\n";
                    ATTR_WARN_J(msg, name, j, uia[j], uib[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_FLOAT) {
            float *fa = (float*) v1->xvalue;
            float *fb = (float*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                /* floating-point inequality here but we genuinely do
                 * expect all processors to set bit-for-bit identical
                 * headers
                if (fa[j] != fb[j]) {
                 */
                if (memcmp(fa+j, fb+j, SIZEOF_FLOAT)) {
                    msg = "%s attribute \"%s\"[%d] FLOAT (root=%f, local=%f)\n";
                    ATTR_WARN_J(msg, name, j, fa[j], fb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_DOUBLE) {
            double *da = (double*) v1->xvalue;
            double *db = (double*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                /* floating-point inequality here but we genuinely do
                 * expect all processors to set bit-for-bit identical
                 * headers
                if (da[j] != db[j]) {
                 */
                if (memcmp(da+j, db+j, SIZEOF_DOUBLE)) {
                    msg = "%s attribute \"%s\"[%d] DOUBLE (root=%f, local=%f)\n";
                    ATTR_WARN_J(msg, name, j, da[j], db[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_INT64) {
            long long *slla = (long long*) v1->xvalue;
            long long *sllb = (long long*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (slla[j] != sllb[j]) {
                    msg = "%s attribute \"%s\"[%d] INT64 (root=%lldll, local=%lldll)\n";
                    ATTR_WARN_J(msg, name, j, slla[j], sllb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        else if (v1->type == NC_UINT64) {
            unsigned long long *ulla = (unsigned long long*) v1->xvalue;
            unsigned long long *ullb = (unsigned long long*) v2->xvalue;
            for (j=0; j<v1->nelems; j++) {
                if (ulla[j] != ullb[j]) {
                    msg = "%s attribute \"%s\"[%d] UINT64 (root=%llull, local=%llull)\n";
                    ATTR_WARN_J(msg, name, j, ulla[j], ullb[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
                    break;
                }
            }
        }
        if (status == NC_NOERR) status = err;

        /* overwrite local's attr with root's */
        if (ErrIsHeaderDiff(err)) {
            ncmpio_free_NC_attr(local_attr->value[i]);
            local_attr->value[i] = ncmpio_dup_NC_attr(root_attr->value[i]);
        }
    }

    /* delete extra attributes defined only in local copy */
    for (; i<local_attr->ndefined; i++)
        ncmpio_free_NC_attr(local_attr->value[i]);

    local_attr->ndefined = root_attr->ndefined;

    return status;
}

/*----< compare_vars() >------------------------------------------------------*/
/* compare the local copy of var_list against root's
 * If inconsistency is detected, overwrite local's with root's
 */
static int
compare_vars(int          safe_mode,
             NC_vararray *root_var,
             NC_vararray *local_var)
{
    int i, j, err, status=NC_NOERR;
    char *msg;

    /* check if the numbers of variables are the same */
    if (root_var->ndefined != local_var->ndefined) {
        if (safe_mode)
            printf("%s number of variables (root=%d, local=%d)\n",
                   WARN_STR, root_var->ndefined, local_var->ndefined);
        DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_VAR_NUM)
    }

    for (i=0; i<root_var->ndefined; i++) {

        if (i >= local_var->ndefined) { /* if local list is shorter */
            /* copy root's variable to local */
            NC_var *new_var = ncmpio_dup_NC_var(root_var->value[i]);
            err = ncmpio_incr_NC_vararray(local_var, new_var);
            if (err != NC_NOERR) return err; /* a fatal error */
            /* local_var->ndefined is increased by 1 in ncmpio_incr_NC_vararray() */
            continue;
        }

        NC_var *v1 = root_var->value[i];
        NC_var *v2 = local_var->value[i];
        char *name = v1->name->cp;

#define VAR_WARN(msg, var, root, local) \
    if (safe_mode) printf(msg, WARN_STR, var, root, local);

#define VAR_WARN_J(msg, var, j, root, local) \
    if (safe_mode) printf(msg, WARN_STR, var, j, root, local);

        err = NC_NOERR;
        if (v1->name->nchars != v2->name->nchars ||
            strncmp(v1->name->cp, v2->name->cp, (size_t)v1->name->nchars) != 0) {
            msg = "%s variable %s (root=%s, local=%s)\n";
            VAR_WARN(msg, "name", name, v2->name->cp)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_NAME)
        }
        else if (v1->ndims != v2->ndims) {
            msg = "%s variable %s's ndims (root=%d, local=%d)\n";
            VAR_WARN(msg, name, v1->ndims, v2->ndims)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_NDIMS)
        }
        else if (v1->type != v2->type) {
            msg = "%s variable %s's type (root=%d, local=%d)\n";
            VAR_WARN(msg, name, v1->type, v2->type)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_TYPE)
        }
        else if (v1->len != v2->len) {
            msg = "%s variable %s's len (root=%lld, local=%lld)\n";
            VAR_WARN(msg, name, v1->len, v2->len)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_LEN)
        }
        else if (v1->begin != v2->begin) {
            msg = "%s variable %s's begin (root=%lld, local=%lld)\n";
            VAR_WARN(msg, name, v1->begin, v2->begin)
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_BEGIN)
        }
        else {
            for (j=0; j<v1->ndims; j++) {
                if (v1->dimids[j] != v2->dimids[j]) {
                    msg = "%s variable %s's %dth dim ID (root=%d, local=%ld)\n";
                    VAR_WARN_J(msg, name, j, v1->dimids[j], v2->dimids[j])
                    DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_VAR_DIMIDS)
                    break;
                }
            }
        }
        /* compare variable's attributes if by far no inconsistency is found */
        if (err == NC_NOERR) {
            err = compare_attrs(safe_mode, &(v1->attrs), &(v2->attrs));
            if (err != NC_NOERR && !ErrIsHeaderDiff(err))
                return err; /* a fatal error */
        }

        if (status == NC_NOERR) status = err;

        /* if there is any inconsistency, overwrite local's var with root's */
        if (ErrIsHeaderDiff(err)) {
            ncmpio_free_NC_var(local_var->value[i]);
            local_var->value[i] = ncmpio_dup_NC_var(root_var->value[i]);
            /* note once a new var is created, one must call
             * compute_var_shape() to recalculate the shape */
        }
    }

    /* delete extra variables defined only in local copy */
    for (; i<local_var->ndefined; i++)
        ncmpio_free_NC_var(local_var->value[i]);

    local_var->ndefined = root_var->ndefined;

#ifndef SEARCH_NAME_LINEARLY
    if (status != NC_NOERR) {
        /* vars are not consistent, must rebuild var name lookup table */
        for (i=0; i<HASH_TABLE_SIZE; i++) {
            if (local_var->nameT[i].num > 0)
                NCI_Free(local_var->nameT[i].list);
            local_var->nameT[i].num = 0;
            local_var->nameT[i].list = NULL;
        }

        /* populate var name lookup table */
        for (i=0; i<local_var->ndefined; i++) {
            /* hash the var name into a key for name lookup */
            int key = HASH_FUNC(local_var->value[i]->name->cp);
            NC_nametable *nameT = &local_var->nameT[key];
            if (nameT->num % NC_NAME_TABLE_CHUNK == 0)
                nameT->list = (int*) NCI_Realloc(nameT->list,
                              (size_t)(nameT->num+NC_NAME_TABLE_CHUNK) * SIZEOF_INT);
            nameT->list[nameT->num] = i;
            nameT->num++;
        }
    }
#endif

    return status;
}

/*----< compare_NC() >--------------------------------------------------------*/
/* This function is only called by NC_check_header()
 * It checks the header of local copy against root's and overwrites local's
 * header object, ncp, with root's header if any inconsistency is detected.
 * This function is called independently and should not contain any MPI
 * communication calls.
 *
 * Possible error codes returned:
 * NC_ENOTNC3, NC_ENOTNC, NC_ESMALL, and all inconsistency errors
 * NC_EMULTIDEFINE_XXX
 */
static int
compare_NC(bufferinfo *getbuf, /* header from root */
           NC         *ncp)
{
    int err, status=NC_NOERR;
    char magic[sizeof(ncmagic1)]; /* root's file format signature */
    MPI_Offset nrecs=0;
    MPI_Aint pos_addr, base_addr;
    NC *root_ncp;

    assert(ncp != NULL);
    assert(getbuf != NULL);

    /* check the file format signature in root's header */
    memset(magic, 0, sizeof(magic));
    err = ncmpix_getn_text((const void **)(&getbuf->pos), sizeof(magic), magic);
    if (err != NC_NOERR) {
        /* Fatal error, as root's header is significant */
        if (ncp->safe_mode) fprintf(stderr,"Error: CDF magic number from root's header\n");
        return err;
    }

    /* check if the first 3 letters are "CDF" */
    if (memcmp(magic, ncmagic1, sizeof(ncmagic1)-1) != 0) {
        /* Fatal error, as root's header is significant */
        /* check if is HDF5 file */
        char signature[8], *hdf5_signature="\211HDF\r\n\032\n";
        memcpy(signature, magic, 4);
        ncmpix_getn_text((const void **)(&getbuf->pos), 4, signature+4);
        if (memcmp(signature, hdf5_signature, 8) == 0) {
            DEBUG_ASSIGN_ERROR(status, NC_ENOTNC3)
            if (ncp->safe_mode)
                fprintf(stderr,"Error: root's header indicates an HDF5 file\n");
        }
        else {
            DEBUG_ASSIGN_ERROR(status, NC_ENOTNC)
            if (ncp->safe_mode)
                fprintf(stderr,"Error: root's header indicates not a CDF file\n");
        }
        return status; /* should not continue */
    }

    /* allocate a header object and fill it with root's header */
    root_ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (root_ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)
    root_ncp->comm    = MPI_COMM_NULL;
    root_ncp->mpiinfo = MPI_INFO_NULL;

    /* in safe_mode, consistency of magic numbers have already been checked in
     * ncmpi_create()
     */
    if (magic[sizeof(ncmagic1)-1] == 0x5) {
        root_ncp->format = 5;
        getbuf->version = 5;
    }
    else if (magic[sizeof(ncmagic1)-1] == 0x2) {
        root_ncp->format = 2;
        getbuf->version = 2;
    }
    else if (magic[sizeof(ncmagic1)-1] != 0x1) {
        getbuf->version = 1;
        root_ncp->format = 1;
        /* Fatal error, as root's header is significant */
        if (ncp->safe_mode)
            fprintf(stderr,"Error: root's header indicates not CDF 1/2/5 format\n");
        DEBUG_RETURN_ERROR(NC_ENOTNC) /* should not continue */
    }

    if (! ncp->safe_mode) {
        /* check local's version number in last byte of magic against root's */
        int root_ver = magic[sizeof(ncmagic1)-1];

        if (ncp->format != root_ver) {
#ifdef PNETCDF_DEBUG
            printf("%s CDF file format (local=CDF-%d, root=CDF-%d)\n",
                   WARN_STR, ncp->format, root_ver);
#endif
            /* overwrite the local header object with root's */
            ncp->format = root_ver;

            /* this inconsistency is not fatal */
            DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_CMODE)
        }
        getbuf->version = root_ver;
        ncp->format     = root_ver;
    }

#if SIZEOF_MPI_OFFSET < 8
    if (getbuf->version > 1) {
        /* for NC_64BIT_DATA or NC_64BIT_OFFSET, MPI_Offset must be 8 bytes */
        if (ncp->safe_mode)
            fprintf(stderr,"Error: cannot support CDF-2 and CDF-5 on this machine\n");
        DEBUG_RETURN_ERROR(NC_ESMALL) /* should not continue */
    }
#endif

    /* since getbuf contains the entire root's header, we do not need to check
     * for any next element in the buffer. Similarly, for all possible calls
     * to hdr_check_buffer() from this subroutine, hdr_fetch() will never be
     * called.
     * (move on to the next element in header: number of records)
     err = hdr_check_buffer(getbuf, (getbuf->version < 5) ? 4 : 8);
     if (err != NC_NOERR) {
         if (ncp->safe_mode)
             fprintf(stderr,"Error: root's header is too short\n");
         return err;
     }
     */

    if (getbuf->version < 5) {
        uint tmp=0;
        err = ncmpix_get_uint32((const void **)(&getbuf->pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    else {
        uint64 tmp=0;
        err = ncmpix_get_uint64((const void **)(&getbuf->pos), &tmp);
        nrecs = (MPI_Offset)tmp;
    }
    if (err != NC_NOERR) {
        if (ncp->safe_mode)
            fprintf(stderr,"Error: failed to read numrecs from root's header\n");
        return err; /* should not continue */
    }

    root_ncp->numrecs = nrecs;
    if (root_ncp->numrecs != ncp->numrecs) {
        if (ncp->safe_mode)
            printf("%s number of records (local=%lld, root=%lld)\n",
                   WARN_STR, ncp->numrecs, root_ncp->numrecs);
        /* overwrite the local header's numrecs */
        ncp->numrecs = root_ncp->numrecs;
        if (status == NC_NOERR) DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE_NUMRECS)
    }

#ifdef HAVE_MPI_GET_ADDRESS
    MPI_Get_address(getbuf->pos,  &pos_addr);
    MPI_Get_address(getbuf->base, &base_addr);
#else
    MPI_Address(getbuf->pos,  &pos_addr);
    MPI_Address(getbuf->base, &base_addr);
#endif
    assert(pos_addr < base_addr + getbuf->size);

    /* get the next header element dim_list from getbuf to root_ncp */
    err = hdr_get_NC_dimarray(getbuf, &root_ncp->dims);
    if (err != NC_NOERR) return err; /* fatal error */

    /* compare local's and root's dim_list */
    err = compare_dims(ncp->safe_mode, &root_ncp->dims, &ncp->dims);
    if (err != NC_NOERR && !ErrIsHeaderDiff(err))
        return err; /* a fatal error */
    if (status == NC_NOERR) status = err;

    /* get the next header element gatt_list from getbuf to root_ncp */
    err = hdr_get_NC_attrarray(getbuf, &root_ncp->attrs);
    if (err != NC_NOERR) return err; /* fatal error */

    /* get the next header element att_list from getbuf to root_ncp */
    err = compare_attrs(ncp->safe_mode, &root_ncp->attrs, &ncp->attrs);
    if (err != NC_NOERR && !ErrIsHeaderDiff(err))
        return err; /* a fatal error */
    if (status == NC_NOERR) status = err;

    /* get the next header element var_list from getbuf to root_ncp */
    err = hdr_get_NC_vararray(getbuf, &root_ncp->vars);
    if (err != NC_NOERR) return err; /* fatal error */

    /* compare local's and root's var_list */
    err = compare_vars(ncp->safe_mode, &root_ncp->vars, &ncp->vars);
    if (err != NC_NOERR && !ErrIsHeaderDiff(err))
        return err; /* a fatal error */
    if (status == NC_NOERR) status = err;

    if (err != NC_NOERR) { /* header has been sync-ed with root */
        /* recompute shape is required for every new variable created */
        err = compute_var_shape(ncp);
        if (err != NC_NOERR) return err; /* a fatal error */
    }
    /* now, the local header object has been sync-ed with root */

    if (ncp->safe_mode && ErrIsHeaderDiff(status)) {
        /* recompute header size */
        MPI_Offset root_xsz, local_xsz;
         root_xsz = ncmpio_hdr_len_NC(root_ncp);
        local_xsz = ncmpio_hdr_len_NC(ncp);
        /* root's header size is getbuf->size */
        assert( ncp->xsz == getbuf->size &&
                root_xsz == local_xsz    &&
               local_xsz == getbuf->size);
    }
    ncmpio_free_NC(root_ncp);
    return status;
}

/*----< NC_check_header() >--------------------------------------------------*/
/*
 * Check the consistency of defined header metadata across all processes and
 * overwrite the local header objects with root's if inconsistency is found.
 * This function is collective.
 */
static int
NC_check_header(NC         *ncp,
                void       *buf,
                MPI_Offset  local_xsz) /* size of buf */
{
    int h_size, rank, g_status, status=NC_NOERR, mpireturn;

    /* root's header size has been broadcasted in NC_begin() and saved in
     * ncp->xsz.
     */

    /* TODO: When root process 0 broadcasts its header,
     * currently the header size cannot be larger than 2^31 bytes,
     * due to the 2nd argument, count, of MPI_Bcast being of type int.
     * Possible solution is to broadcast in chunks of 2^31 bytes.
     */
    h_size = (int)ncp->xsz;
    if (ncp->xsz != h_size)
        DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

    MPI_Comm_rank(ncp->nciop->comm, &rank);

    if (rank == 0) {
        TRACE_COMM(MPI_Bcast)(buf, h_size, MPI_BYTE, 0, ncp->nciop->comm);
    }
    else {
        bufferinfo gbp;
        void *cmpbuf = (void*) NCI_Malloc((size_t)h_size);

        TRACE_COMM(MPI_Bcast)(cmpbuf, h_size, MPI_BYTE, 0, ncp->nciop->comm);

        if (h_size != local_xsz || memcmp(buf, cmpbuf, h_size)) {
            /* now part of this process's header is not consistent with root's
             * check and report the inconsistent part
             */

            /* Note that gbp.nciop and gbp.offset below will not be used in
             * compare_NC() */
            gbp.nciop  = ncp->nciop;
            gbp.offset = 0;
            gbp.size   = h_size;   /* entire header is in the buffer, cmpbuf */
            gbp.index  = 0;
            gbp.pos    = gbp.base = cmpbuf;

            /* find the inconsistent part of the header, report the difference,
             * and overwrite the local header object with root's.
             * compare_NC() should not have any MPI communication calls.
             */
            status = compare_NC(&gbp, ncp);

            /* header consistency is only checked on non-root processes. The
             * returned status can be a fatal error or header inconsistency
             * error, (fatal errors are due to object allocation), but never
             * NC_NOERR.
             */
        }
        NCI_Free(cmpbuf);
    }

    if (ncp->safe_mode) {
        TRACE_COMM(MPI_Allreduce)(&status, &g_status, 1, MPI_INT, MPI_MIN,
                                  ncp->nciop->comm);
        if (mpireturn != MPI_SUCCESS) {
            return ncmpio_handle_error(mpireturn, "MPI_Allreduce");
        }
        if (g_status != NC_NOERR) { /* some headers are inconsistent */
            if (status == NC_NOERR) DEBUG_ASSIGN_ERROR(status, NC_EMULTIDEFINE)
        }
    }

    return status;
}
#endif

/*----< ncmpio_read_NC() >---------------------------------------------------*/
/* Re-read in the header.
 * On PnetCDF, this is of no use, as header metadata is always sync-ed among
 * all processes, except for numrecs, which can be sync-ed by calling
 * ncmpio_sync_numrecs()
 */
int
ncmpio_read_NC(NC *ncp) {
    int status = NC_NOERR;

    ncmpio_free_NC_dimarray(&ncp->dims);
    ncmpio_free_NC_attrarray(&ncp->attrs);
    ncmpio_free_NC_vararray(&ncp->vars);

    status = ncmpio_hdr_get_NC(ncp);

    if (status == NC_NOERR)
        fClr(ncp->flags, NC_NDIRTY);

    return status;
}

/*----< ncmpio_find_NC_Udim() >----------------------------------------------*/
/*
 * Step thru NC_DIMENSION array, seeking the UNLIMITED dimension.
 * Return dimid or -1 on not found.
 * *dimpp is set to the appropriate NC_dim.
 */
int
ncmpio_find_NC_Udim(const NC_dimarray  *ncap,
                    NC_dim            **dimpp)
{
    int dimid;

    assert(ncap != NULL);

    if (ncap->ndefined == 0) return -1;

    /* note that the number of dimensions allowed is < 2^32 */
    for (dimid=0; dimid<ncap->ndefined; dimid++)
        if (ncap->value[dimid]->size == NC_UNLIMITED) {
            /* found the matched name */
            if (dimpp != NULL)
                *dimpp = ncap->value[dimid];
            return dimid;
        }

    /* not found */
    return -1;
}

/*----< ncmpio_swapn() >-----------------------------------------------------*/
/* out-place byte swap, i.e. dest_p != src_p */
void
ncmpio_swapn(void       *dest_p,  /* destination array */
             const void *src_p,   /* source array */
             MPI_Offset  nelems,  /* number of elements in buf[] */
             int         esize)   /* byte size of each element */
{
    int  i;

    if (esize <= 1 || nelems <= 0) return;  /* no need */

    if (esize == 4) { /* this is the most common case */
              uint32_t *dest = (uint32_t*)       dest_p;
        const uint32_t *src  = (const uint32_t*) src_p;
        for (i=0; i<nelems; i++) {
            dest[i] = src[i];
            dest[i] =  ((dest[i]) << 24)
                    | (((dest[i]) & 0x0000ff00) << 8)
                    | (((dest[i]) & 0x00ff0000) >> 8)
                    | (((dest[i]) >> 24));
        }
    }
    else if (esize == 8) {
              uint64_t *dest = (uint64_t*)       dest_p;
        const uint64_t *src  = (const uint64_t*) src_p;
        for (i=0; i<nelems; i++) {
            dest[i] = src[i];
            dest[i] = ((dest[i] & 0x00000000000000FFULL) << 56) |
                      ((dest[i] & 0x000000000000FF00ULL) << 40) |
                      ((dest[i] & 0x0000000000FF0000ULL) << 24) |
                      ((dest[i] & 0x00000000FF000000ULL) <<  8) |
                      ((dest[i] & 0x000000FF00000000ULL) >>  8) |
                      ((dest[i] & 0x0000FF0000000000ULL) >> 24) |
                      ((dest[i] & 0x00FF000000000000ULL) >> 40) |
                      ((dest[i] & 0xFF00000000000000ULL) >> 56);
        }
    }
    else if (esize == 2) {
              uint16_t *dest =       (uint16_t*) dest_p;
        const uint16_t *src  = (const uint16_t*) src_p;
        for (i=0; i<nelems; i++) {
            dest[i] = src[i];
            dest[i] = (uint16_t)(((dest[i] & 0xff) << 8) |
                                 ((dest[i] >> 8) & 0xff));
        }
    }
    else {
              uchar *op = (uchar*) dest_p;
        const uchar *ip = (uchar*) src_p;
        /* for esize is not 1, 2, or 4 */
        while (nelems-- > 0) {
            for (i=0; i<esize; i++)
                op[i] = ip[esize-1-i];
            op += esize;
            ip += esize;
        }
    }
}

inline nc_type
ncmpio_mpi2nctype(MPI_Datatype itype)
{
    if (itype == MPI_SIGNED_CHAR)        return NC_BYTE ;
    if (itype == MPI_CHAR)               return NC_CHAR ;
    if (itype == MPI_SHORT)              return NC_SHORT ;
    if (itype == MPI_INT)                return NC_INT ;
    if (itype == MPI_FLOAT)              return NC_FLOAT ;
    if (itype == MPI_DOUBLE)             return NC_DOUBLE ;
    if (itype == MPI_UNSIGNED_CHAR)      return NC_UBYTE ;
    if (itype == MPI_UNSIGNED_SHORT)     return NC_USHORT ;
    if (itype == MPI_UNSIGNED)           return NC_UINT ;
    if (itype == MPI_LONG_LONG_INT)      return NC_INT64 ;
    if (itype == MPI_UNSIGNED_LONG_LONG) return NC_UINT64 ;
    return NC_EBADTYPE;
}

/*----< ncmpi_print_all_var_offsets() >---------------------------------------*/
/* This is an independent subroutine */
int
ncmpi_print_all_var_offsets(int ncid) {
    int i, err;
    NC_var **vpp=NULL;
    NC *ncp=NULL;
    PNC *pncp;

    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
    ncp = (NC*)pncp->ncp;

    if (ncp->begin_var%1048576)
        printf("%s header size (ncp->begin_var)=%lld MB + %lld\n",
        ncp->path, ncp->begin_var/1048575, ncp->begin_var%1048576);
    else
        printf("%s header size (ncp->begin_var)=%lld MB\n",
        ncp->path, ncp->begin_var/1048575);

    vpp = ncp->vars.value;
    for (i=0; i<ncp->vars.ndefined; i++, vpp++) {
        char str[1024];
        MPI_Offset off = (*vpp)->begin;
        MPI_Offset rem = off % 1048576;;

        if (IS_RECVAR(*vpp))
            sprintf(str,"    Record variable \"%20s\": ",(*vpp)->name->cp);
        else
            sprintf(str,"non-record variable \"%20s\": ",(*vpp)->name->cp);

        if (rem)
            printf("%s offset=%12lld MB + %7lld len=%lld\n", str, off/1048576, rem,(*vpp)->len);
        else
            printf("%s offset=%12lld MB len=%lld\n", str, off/1048576,(*vpp)->len);
    }
    return NC_NOERR;
}

/*----< ncmpio_free_NC_string() >--------------------------------------------*/
/*
 * Free string, and, if needed, its values.
 * Formerly
NC_free_string()
 */
inline void
ncmpio_free_NC_string(NC_string *ncstrp)
{
    if (ncstrp == NULL) return;
    /* ncstrp->cp is allocated as part of ncstrp object */
    NCI_Free(ncstrp);
}

/*----< ncmpio_new_NC_string() >---------------------------------------------*/
/*
 * Allocate a NC_string structure large enough
 * to hold slen characters.
 */
NC_string *
ncmpio_new_NC_string(size_t      slen,
                     const char *str)  /* must be already normalized */
{
    /* str may not be NULL terminated */
    NC_string *ncstrp;
    size_t sizeof_NC_string = M_RNDUP(sizeof(NC_string));
    size_t sz = slen + sizeof_NC_string + 1;
    /* one char more space for NULL terminate char */

    ncstrp = (NC_string *) NCI_Calloc(sz, sizeof(char));
    if (ncstrp == NULL) return NULL;

    /* make space occupied by ncstrp->cp part of ncstrp */
    ncstrp->nchars = (MPI_Offset)slen;
    ncstrp->cp = (char *)ncstrp + sizeof_NC_string;

    /* in PnetCDF, we want to make name->cp always NULL character terminated */
    if (str != NULL && *str != '\0') {
        strncpy(ncstrp->cp, str, slen);
        ncstrp->cp[slen] = '\0';  /* NULL terminated */
    }

    return(ncstrp);
}


/*----< ncmpio_set_NC_string() >---------------------------------------------*/
/*
 * If possible, change the value of an NC_string to 'str'.
 */
int
ncmpio_set_NC_string(NC_string  *ncstrp,
                     const char *str)
{
    size_t slen;

    assert(str != NULL && *str != '\0');

    slen = strlen(str);

    if (ncstrp->nchars < (MPI_Offset)slen) DEBUG_RETURN_ERROR(NC_ENOTINDEFINE)

    memcpy(ncstrp->cp, str, slen);

    /* in PnetCDF, we want to make name->cp always NULL character terminated */
    if (ncstrp->nchars > (MPI_Offset)slen)
        memset(ncstrp->cp + slen, 0, (size_t)ncstrp->nchars - slen);

    ncstrp->nchars = (MPI_Offset)slen;

    return NC_NOERR;
}

/*@
  ncmpio_data_repack - copy data between two buffers with different datatypes.

  Input:
. inbuf - input buffer where data is copied from
. incount - number of input elements
. intype - datatype of each element in input buffer
. outbuf - output buffer where data is copied to
. outcount - number of output elements
. outtype - datatype of each element in output buffer
@*/

int ncmpio_data_repack(void *inbuf,
                       MPI_Offset incount,
                       MPI_Datatype intype,
                       void *outbuf,
                       MPI_Offset outcount,
                       MPI_Datatype outtype)
{
  int intypesz, outtypesz;
  int packsz;
  void *packbuf;
  int packpos;

  MPI_Type_size(intype, &intypesz);
  MPI_Type_size(outtype, &outtypesz);

  if (incount*intypesz != outcount*outtypesz) {
    /* input data amount does not match output data amount */
    /* NOTE: we ignore it for user responsibility or add error handling ? */

    /* for rescue, guarantee output data amount <= input data amount */
    if (incount*intypesz < outcount*outtypesz)
      outcount = incount*intypesz/outtypesz;
  }

  if (incount == 0)
    return NC_NOERR;

  /* local pack-n-unpack, using MPI_COMM_SELF */
  if (incount  != (int)incount)  DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
  if (outcount != (int)outcount) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)

  MPI_Pack_size((int)incount, intype, MPI_COMM_SELF, &packsz);
  packbuf = (void *)NCI_Malloc((size_t)packsz);
  packpos = 0;
  MPI_Pack(inbuf, (int)incount, intype, packbuf, packsz, &packpos, MPI_COMM_SELF);
  packpos = 0;
  MPI_Unpack(packbuf, packsz, &packpos, outbuf, (int)outcount, outtype, MPI_COMM_SELF);
  NCI_Free(packbuf);

  return NC_NOERR;
}

/*----< ncmpio_cktype() >----------------------------------------------------*/
/*  Verify that this is a user nc_type */
int
ncmpio_cktype(int cdf_ver, nc_type type)
{
    /* the max data type supported by CDF-5 is NC_UINT64 */
    if (type <= 0 || type > NC_UINT64)
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    /* For CDF-1 and CDF-2 files, only classic types are allowed. */
    if (cdf_ver < 5 && type > NC_DOUBLE)
        DEBUG_RETURN_ERROR(NC_ESTRICTCDF2)

    return NC_NOERR;
}

