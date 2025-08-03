#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>
#include "baseline_ncx_app.h" 
#include "mem_tracker.h"

#ifdef MEM_TRACKING
#define malloc(size) tracked_malloc(size)
#define free(ptr)    tracked_free(ptr)
#define realloc(ptr, size) tracked_realloc(ptr, size)
#define strdup(s) tracked_strdup(s)
#endif


  /* ---------------------------------- Serializaition ----------------------------------------*/

int
xlen_nc_type_meta(nc_type xtype, int *size)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  *size = 1; return 0;
        case NC_SHORT:
        case NC_USHORT: *size = 2; return 0;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  *size = 4; return 0;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: *size = 8; return 0;
    }
    return 0;
}

static int 
putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
	(void) memcpy(*xpp, tp, (size_t)nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return 0;
}

static int
put_uint32(void **xpp, unsigned int ip)
{
    memcpy(*xpp, &ip, 4);
    /* advace *xpp 4 bytes */
    *xpp  = (void *)((char *)(*xpp) + 4);
    return 0;
}

static int
serialize_dim(metabuffer   *pbp,
               const hdr_dim *dimp)
{
    /* copy name */
    serialize_name(pbp, dimp->name);
    put_uint32((void**)(&pbp->pos), (uint32_t)dimp->size);
    return 0;
}

static int
serialize_name(metabuffer *pbp,
                const char *name)
{
    size_t nchars = strlen(name);

    put_uint32((void**)(&pbp->pos), (uint32_t)nchars);

    return putn_text((void **)(&pbp->pos), (MPI_Offset)nchars, name);
}


static int
serialize_dimarray(metabuffer        *pbp,
                    const hdr_dimarray *ncap)
{
    int i, status;
    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_DIMENSION */
        status = put_uint32((void**)(&pbp->pos), NC_DIMENSION);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;
        /* copy [dim ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = serialize_dim(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }
    return 0;
}

static int
serialize_attrV(metabuffer    *pbp,
                 const hdr_attr *attrp)
{

    int xsz;
    int sz;

    /* xlen_nc_type_meta() returns the element size (unaligned) of
     * attrp->xtype attrp->xsz is the aligned total size of attribute values
     */
    xlen_nc_type_meta(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;
    memcpy(pbp->pos, attrp->xvalue, (size_t)sz);
    pbp->pos = (void *)((char *)pbp->pos + sz);
    return 0;
}

/*----< serialize_NC_attr() >--------------------------------------------------*/
static int
serialize_attr(metabuffer    *pbp,
                const hdr_attr *attrp)
{
    int status;

    /* copy name */
    status = serialize_name(pbp, attrp->name);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)attrp->xtype);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)attrp->nelems);
    if (status != NC_NOERR) return status;

    /* copy [values ...] */
    status = serialize_attrV(pbp, attrp);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

static int
serialize_attrarray(metabuffer         *pbp,
                     const hdr_attrarray *ncap)
{

    int i, status;
    assert(pbp != NULL);
    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_ATTRIBUTE */
        status = put_uint32((void**)(&pbp->pos), NC_ATTRIBUTE);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;
        /* copy [attr ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = serialize_attr(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

static int
serialize_var(metabuffer   *pbp,
               const hdr_var *varp)
{
    int i, status;

    /* copy name */
    status = serialize_name(pbp, varp->name);
    if (status != NC_NOERR) return status;

    /* copy nelems */

    status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->ndims);

    if (status != NC_NOERR) return status;

    /* copy [dim_index ...] i*/
    for (i=0; i<varp->ndims; i++) {

        status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->dimids[i]);
        if (status != NC_NOERR) return status;
    }

    /* copy vatt_list */
    status = serialize_attrarray(pbp, &varp->attrs);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = put_uint32((void**)(&pbp->pos), (uint32_t)varp->xtype);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}


/*----< serialize_vararray() >----------------------------------------------*/
static int
serialize_vararray(metabuffer        *pbp,
                    const hdr_vararray *ncap)
{
    int i, status;
    assert(pbp != NULL);
    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;
        status = put_uint32((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_VARIABLE */
        status = put_uint32((void**)(&pbp->pos), NC_VARIABLE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = put_uint32((void**)(&pbp->pos), (uint32_t)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [var ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status =serialize_var(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }
    return NC_NOERR;
}



/*----< serialize_hdr_meta() >----------------------------------------------*/
int
serialize_hdr_meta(struct hdr *ncp, void *buf)
{
    int status;
    metabuffer putbuf;

    putbuf.pos           = buf;
    putbuf.base          = buf;
    putbuf.size          = ncp->xsz;

    /* copy dim_list */
    status = serialize_dimarray(&putbuf, &ncp->dims);
    if (status != NC_NOERR) return status;



    // /* copy gatt_list */
    // status = serialize_attrarray(&putbuf, &ncp->attrs);
    // if (status != NC_NOERR) return status;

    /* copy var_list */
    status = serialize_vararray(&putbuf, &ncp->vars);
    if (status != NC_NOERR) return status;

    size_t serializedSize = putbuf.pos - putbuf.base;

    // Print the result
    // printf("Number of bytes taken after serialization: %zu\n", serializedSize);



    return NC_NOERR;
}

  /* ---------------------------------- Deserializaition ----------------------------------------*/


static int
getn_text(void **xpp, MPI_Offset nelems, char *tp)
{
	(void) memcpy(tp, *xpp, (size_t)nelems);
    tp[nelems] = '\0';
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}


static int
get_uint32(void **xpp, unsigned int *ip)
{
    memcpy(ip, *xpp, 4);
    /* advance *xpp 4 bytes */
    *xpp = (void *)((const char *)(*xpp) + 4);
    return NC_NOERR;
}

static int deserialize_nc_type(metabuffer *gbp, nc_type *xtypep){
    int err;
    uint32_t xtype;
    err = get_uint32((void **)(&gbp->pos), &xtype);
    if (err != NC_NOERR) return err;
    *xtypep = (nc_type) xtype;
    return NC_NOERR;
}

static int deserialize_name(metabuffer *gbp, char **name) {
    unsigned int nchars;
    get_uint32((void**)&gbp->pos, &nchars);
    *name = (char *)malloc(nchars + 1);
    if (*name == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }
    getn_text((void **)&gbp->pos, nchars, *name);
    return NC_NOERR;
}

static int deserialize_dim(metabuffer *gbp, hdr_dim *dimp) {
    MPI_Offset dim_length;
    uint32_t tmp;
    char *name;
    int err;
    err = deserialize_name(gbp, &name); 
    if (err != NC_NOERR) return err;
    get_uint32((void**)&gbp->pos, &tmp);
    dim_length = (MPI_Offset)tmp;
    dimp->name     = name;
    dimp->name_len = strlen(name);
    dimp->size     = dim_length;
    return 0;
}

static int deserialize_dimarray(metabuffer *gbp, hdr_dimarray *ncap) {
    unsigned int tag;
    get_uint32((void**)&gbp->pos, &tag);
    if (tag == NC_UNSPECIFIED) {
        get_uint32((void**)&gbp->pos, (unsigned int *)&ncap->ndefined);
        assert(ncap->ndefined == 0);
        return 0; // ABSENT
    } else if (tag == NC_DIMENSION) {
        get_uint32((void**)&gbp->pos, (unsigned int *)&ncap->ndefined);

        ncap->value = (hdr_dim **)malloc(ncap->ndefined * sizeof(hdr_dim *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }

        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_dim *)malloc(sizeof(hdr_dim));
            if (ncap->value[i] == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                return -1;
            }
            if (deserialize_dim(gbp, ncap->value[i]) != 0) {
                return -1;
            }
        }
    }

    return 0;
}

static int deserialize_attrV(metabuffer *gbp, hdr_attr *attrp) {
    int xsz, sz, err;

    xlen_nc_type_meta(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;

    attrp->xvalue = malloc(sz);
    if (attrp->xvalue == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }

    memcpy(attrp->xvalue, gbp->pos, (size_t)sz);
    gbp->pos = (void *)((char *)(gbp->pos) + sz);

    return 0;
}

static int deserialize_attr(metabuffer *gbp, hdr_attr *attrp) {
    uint32_t tmp;
    int err;
    char *name;
    err = deserialize_name(gbp, &name);
    if (err != NC_NOERR) return err;
    attrp->name = name;
    attrp->name_len = strlen(name);
    err = deserialize_nc_type(gbp, &attrp->xtype);
    if (err != NC_NOERR) return err;
    err = get_uint32((void**)&gbp->pos, &tmp);
    attrp->nelems = (int)tmp;
    if (err != NC_NOERR) return err;
    err = deserialize_attrV(gbp, attrp);
    if (err != NC_NOERR) return err;

    return 0;
}

static int deserialize_attrarray(metabuffer *gbp, hdr_attrarray *ncap) {
    unsigned int tag;
    get_uint32((void**)&gbp->pos, &tag);
    uint32_t tmp;

    if (tag == NC_UNSPECIFIED) {
        get_uint32((void**)&gbp->pos, &tmp);
        ncap->ndefined = (int) tmp;
        assert(ncap->ndefined == 0);
        return 0; // ABSENT
    } else if (tag == NC_ATTRIBUTE) {
        get_uint32((void**)&gbp->pos, &tmp);
        ncap->ndefined = (int) tmp;
        ncap->value = (hdr_attr **)malloc(ncap->ndefined * sizeof(hdr_attr *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }
        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_attr *)malloc(sizeof(hdr_attr));
            if (ncap->value[i] == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                return -1;
            }

            if (deserialize_attr(gbp, ncap->value[i]) != 0) {
                return -1;
            }
        }
    }

    return 0;
}

static int deserialize_var(metabuffer *gbp, hdr_var *varp) {
    int err;
    char *name;
    // if (deserialize_name(gbp, &varp->name) != 0) {
    //     return -1;
    // }
    /* get name */
    varp->name = NULL;
    varp->name_len = 0;
    varp->dimids = NULL;
    varp->attrs.ndefined = 0;
    varp->attrs.value = NULL;
    err = deserialize_name(gbp, &name);
    if (err != NC_NOERR) return err;
    varp->name = name;
    varp->name_len = strlen(name);
    /* nelems (number of dimensions) */
    u_int32_t tmp;
    get_uint32((void**)&gbp->pos, (unsigned int *)&tmp);
    varp->ndims = (int) tmp;
    varp->dimids = (int *)malloc(varp->ndims * sizeof(int));
    if (varp->dimids == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }

    for (int i = 0; i < varp->ndims; i++) {
        get_uint32((void**)&gbp->pos, &tmp);
        varp->dimids[i] = (int)tmp;
    }

    if (deserialize_attrarray(gbp, &varp->attrs) != 0) {
        return -1;
    }
    err = deserialize_nc_type(gbp, &varp->xtype);
    if (err != NC_NOERR) return err;

    return 0;
}

static int deserialize_vararray(metabuffer *gbp, hdr_vararray *ncap) {
    unsigned int tag;
    int err;
    uint32_t tmp;
    get_uint32((void**)&gbp->pos, &tag);

    if (tag == NC_UNSPECIFIED) {
        get_uint32((void**)&gbp->pos, (unsigned int *)&tmp);
        ncap->ndefined = (int)tmp;
        assert(ncap->ndefined == 0);
        return 0; // ABSENT
    } else if (tag == NC_VARIABLE) {
        get_uint32((void**)&gbp->pos, (unsigned int *)&tmp);
        ncap->ndefined = (int)tmp;
        ncap->value = (hdr_var **)malloc(ncap->ndefined * sizeof(hdr_var *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }

        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_var *)malloc(sizeof(hdr_var));
            if (ncap->value[i] == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                return -1;
            }
            if (deserialize_var(gbp, ncap->value[i]) != 0) {
                return -1;
            }
        }
    }

    return 0;
}

int deserialize_hdr_meta(struct hdr *ncp, void *buf, int buf_size) {

    int status;
    metabuffer getbuf;

    getbuf.pos           = buf;
    getbuf.base          = buf;
    getbuf.size          = buf_size;

    /* get dim_list from getbuf into ncp */
    status = deserialize_dimarray(&getbuf, &ncp->dims);
    if (status != NC_NOERR) return status;
    

    status = deserialize_vararray(&getbuf, &ncp->vars);
    if (status != NC_NOERR) return status;
    // printf("HERE: %ld", getbuf.pos - getbuf.base);
    // printf("HERE: %ld",  getbuf.size);
    // assert((int)(getbuf.pos - getbuf.base) == getbuf.size);
    if ((int)(getbuf.pos - getbuf.base) != getbuf.size) {
        printf("Deserialization error: consumed = %ld, expected = %d\n",
            (long)(getbuf.pos - getbuf.base), getbuf.size);
        return -1;
    }



    return 0;
}



void free_hdr_dim_meta(hdr_dim *dim) {
    if (dim != NULL) {
        free(dim->name);
        free(dim);
    }
}

void free_hdr_dimarray_meta(hdr_dimarray *dims) {
    if (dims != NULL) {
        for (int i = 0; i < dims->ndefined; i++) {
            
            free_hdr_dim_meta(dims->value[i]);
        }
        free(dims->value);
        //free(dims);
    }
}

void free_hdr_attr_meta(hdr_attr *attr) {
    if (attr != NULL) {
        free(attr->name);
        free(attr->xvalue);
    }
}

void free_hdr_attrarray_meta(hdr_attrarray *attrs) {
    if (attrs != NULL) {
        if (attrs->value != NULL) {
            for (int i = 0; i < attrs->ndefined; i++) {
                free_hdr_attr_meta(attrs->value[i]);
            }
            free(attrs->value);
            attrs->value = NULL;
            // free(attrs);
        }
    }
}

void free_hdr_var_meta(hdr_var *var) {
    if (var != NULL) {
        free(var->name);
        free(var->dimids);

        free_hdr_attrarray_meta(&(var->attrs));
        free(var);
    }
}

void free_hdr_vararray_meta(hdr_vararray *vars) {
    if (vars != NULL) {
        for (int i = 0; i < vars->ndefined; i++) {
            free_hdr_var_meta(vars->value[i]);
        }
        free(vars->value);
        // free(vars);
    }
}

void free_hdr_meta(struct hdr *header) {
    if (header != NULL) {
        free_hdr_dimarray_meta(&(header->dims));
        // free_hdr_attrarray_meta(&(header->attrs));
        free_hdr_vararray_meta(&(header->vars));
    }
}

/* ---------------------------------- Read Metadata ----------------------------------------*/
int read_metadata_from_file(const char* filename, struct hdr *recv_hdr) {
    // Open the file in binary read mode
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

    // Move file pointer to the end to determine the file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);  // Move the file pointer back to the beginning

    // Allocate memory for the buffer
    // char *buffer = (char*) malloc(file_size);
    char *buffer = (char*) malloc(file_size);
    if (buffer == NULL) {
        perror("Failed to allocate memory");
        fclose(file);
        return 1;
    }

    // Read the file contents into the buffer
    // printf("file_size: %ld\n", file_size);
    size_t read_size = fread(buffer, 1, file_size, file);
    if (read_size != file_size) {
        perror("Failed to read complete file");
        free(buffer);
        fclose(file);
        return 1;
    }

    // Close the file after reading
    fclose(file);
    recv_hdr->xsz = file_size;
    deserialize_hdr_meta(recv_hdr, buffer, file_size);
    // free(buffer);
    free(buffer);
    return 0;
}
  
/* ---------------------------------- Distribute Metadata ----------------------------------------*/

int distribute_metadata(int rank, int nproc, struct hdr *all_hdr, struct hdr *local_hdr) {
    int num_vars = all_hdr->vars.ndefined;
    int vars_per_process = num_vars / nproc;
    int remainder = num_vars % nproc;

    // Determine start and count based on rank
    int start = rank * vars_per_process + (rank < remainder ? rank : remainder);
    int count = vars_per_process + (rank < remainder ? 1 : 0);

    local_hdr->vars.ndefined = count;
    local_hdr->dims.ndefined = 0;
    local_hdr->attrs.ndefined = 0;
    local_hdr->xsz = 2 * sizeof(uint32_t); // NC_Variable and ndefined
    local_hdr->xsz += 2 * sizeof(uint32_t); // NC_Dimension and nelems

    local_hdr->vars.value = (hdr_var **)malloc(count * sizeof(hdr_var *));
    local_hdr->dims.value = NULL;

    int tot_num_dims = 0;

    for (int i = 0; i < count; ++i) {
        hdr_var *src_var = all_hdr->vars.value[start + i];
        hdr_var *dst_var = (hdr_var *)malloc(sizeof(hdr_var));

        // Copy name
        dst_var->name_len = src_var->name_len;
        dst_var->name = strdup(src_var->name);

        // Copy type and dims
        dst_var->xtype = src_var->xtype;
        dst_var->ndims = src_var->ndims;
        dst_var->dimids = (int *)malloc(src_var->ndims * sizeof(int));

        // Allocate and copy attributes
        dst_var->attrs.ndefined = src_var->attrs.ndefined;
        dst_var->attrs.value = (hdr_attr **)malloc(src_var->attrs.ndefined * sizeof(hdr_attr *));

        for (int j = 0; j < src_var->attrs.ndefined; ++j) {
            hdr_attr *src_attr = src_var->attrs.value[j];
            hdr_attr *dst_attr = (hdr_attr *)malloc(sizeof(hdr_attr));
            dst_attr->name_len = src_attr->name_len;
            dst_attr->name = strdup(src_attr->name);
            dst_attr->xtype = src_attr->xtype;
            dst_attr->nelems = src_attr->nelems;

            int elem_sz;
            xlen_nc_type_meta(dst_attr->xtype, &elem_sz);
            dst_attr->xvalue = malloc(elem_sz * dst_attr->nelems);
            memcpy(dst_attr->xvalue, src_attr->xvalue, elem_sz * dst_attr->nelems);

            dst_var->attrs.value[j] = dst_attr;

            local_hdr->xsz += sizeof(uint32_t) + dst_attr->name_len * sizeof(char);
            local_hdr->xsz += 2 * sizeof(uint32_t);
            local_hdr->xsz += dst_attr->nelems * elem_sz;
        }

        // Copy dimids and related dim info
        for (int k = 0; k < dst_var->ndims; ++k) {
            int dimid = src_var->dimids[k];
            hdr_dim *src_dim = all_hdr->dims.value[dimid];

            // Copy dimension
            hdr_dim *dst_dim = (hdr_dim *)malloc(sizeof(hdr_dim));
            dst_dim->name_len = src_dim->name_len;
            dst_dim->name = strdup(src_dim->name);
            dst_dim->size = src_dim->size;

            // Reallocate dims array
            local_hdr->dims.value = (hdr_dim **)realloc(local_hdr->dims.value, (local_hdr->dims.ndefined + 1) * sizeof(hdr_dim *));
            local_hdr->dims.value[local_hdr->dims.ndefined] = dst_dim;

            dst_var->dimids[k] = local_hdr->dims.ndefined;
            local_hdr->dims.ndefined++;

            local_hdr->xsz += sizeof(uint32_t) + dst_dim->name_len * sizeof(char);
            local_hdr->xsz += sizeof(uint32_t);
        }

        local_hdr->vars.value[i] = dst_var;
        local_hdr->xsz += sizeof(uint32_t) + dst_var->name_len * sizeof(char); // var name
        local_hdr->xsz += sizeof(uint32_t); // xtype
        local_hdr->xsz += sizeof(uint32_t); // nelems of dim list
        local_hdr->xsz += sizeof(uint32_t) * dst_var->ndims; // dimid list
        local_hdr->xsz += 2 * sizeof(uint32_t); // NC_Attribute and ndefined
    }
}