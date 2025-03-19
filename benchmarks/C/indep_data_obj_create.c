/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#include <math.h>


#include <assert.h>



static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define OUTPUT_NAME "benchmark_testfile.nc"
#define NUM_VARS 10
#define VAR_TYPE NC_INT
#define NUM_DIMS 2
#define NUM_ATTRS 0
#define DIM_SIZE 10
#define ATTR_SIZE 10
#define HASH_SIZE 4096

double def_start_time, total_def_time = 0;



#define NC_NOERR        0  



// Static variables to track memory usage
static size_t current_memory_usage = 0;
static size_t max_memory_usage = 0;

// Hash table constants
#define HASH_TABLE_SIZE (1024 * 64)
#define CHUNK_SIZE 1024  // Array grows by this size

// Allocation metadata
typedef struct {
    void* ptr;
    size_t size;
} AllocationNode;

// Bucket for handling collisions
typedef struct {
    AllocationNode* nodes;  // Array of allocation nodes
    size_t capacity;        // Current capacity of the array
    size_t count;           // Number of entries in the array
} Bucket;

// Hash table for allocation tracking
static Bucket hash_table[HASH_TABLE_SIZE] = {0};

// Hash function for pointers
static size_t hash_ptr(void* ptr) {
    uintptr_t address = (uintptr_t)ptr;  // Convert pointer to an integer
    address = (address >> 4) ^ (address << 8);  // Mix bits by shifting
    return (address * 2654435761u) % HASH_TABLE_SIZE;
}

// Add allocation to the hash table
static void add_allocation(void* ptr, size_t size) {
    size_t index = hash_ptr(ptr);
    Bucket* bucket = &hash_table[index];

    // Allocate memory for the bucket if it is empty
    if (bucket->nodes == NULL) {
        bucket->nodes = (AllocationNode*)calloc(CHUNK_SIZE, sizeof(AllocationNode));
        if (!bucket->nodes) {
            fprintf(stderr, "Error: Memory allocation failed for bucket.\n");
            exit(EXIT_FAILURE);
        }
        bucket->capacity = CHUNK_SIZE;
        bucket->count = 0;
    }

    // Resize the bucket array if needed
    if (bucket->count >= bucket->capacity) {
        size_t new_capacity = bucket->capacity + CHUNK_SIZE;
        AllocationNode* new_nodes = (AllocationNode*)realloc(bucket->nodes, new_capacity * sizeof(AllocationNode));
        if (!new_nodes) {
            fprintf(stderr, "Error: Memory reallocation failed for bucket.\n");
            exit(EXIT_FAILURE);
        }
        bucket->nodes = new_nodes;
        bucket->capacity = new_capacity;
    }

    // Add the new allocation to the bucket
    bucket->nodes[bucket->count].ptr = ptr;
    bucket->nodes[bucket->count].size = size;
    bucket->count++;
}

// Remove allocation from the hash table and return its size
static size_t remove_allocation(void* ptr) {
    size_t index = hash_ptr(ptr);
    Bucket* bucket = &hash_table[index];

    if (bucket->nodes == NULL) {
        fprintf(stderr, "Warning: Attempt to free untracked pointer.\n");
        return 0; // Bucket is empty
    }

    // Find the allocation in the bucket
    for (size_t i = 0; i < bucket->count; i++) {
        if (bucket->nodes[i].ptr == ptr) {
            size_t size = bucket->nodes[i].size;

            // Replace the removed entry with the last entry for O(1) removal
            bucket->nodes[i] = bucket->nodes[bucket->count - 1];
            bucket->count--;

            return size;
        }
    }

    fprintf(stderr, "Warning: Attempt to free untracked pointer.\n");
    return 0; // Pointer not found
}

// Free all allocations in the hash table
void free_allocation_struct(void) {
    for (size_t i = 0; i < HASH_TABLE_SIZE; i++) {
        Bucket* bucket = &hash_table[i];
        if (bucket->nodes) {
            free(bucket->nodes); // Free the bucket's dynamic array
            bucket->nodes = NULL;
            bucket->capacity = 0;
            bucket->count = 0;
        }
    }
}

// Wrapper for malloc
void* tracked_malloc(size_t size) {
    void* ptr = malloc(size);
    if (ptr) {
        add_allocation(ptr, size);
        current_memory_usage += size;
        if (current_memory_usage > max_memory_usage) {
            max_memory_usage = current_memory_usage;
        }
    }
    return ptr;
}

// Wrapper for free
void tracked_free(void* ptr) {
    if (ptr) {
        size_t size = remove_allocation(ptr);
        current_memory_usage -= size;
        free(ptr);
    }
}

// Returns current memory usage
size_t inq_malloc_use(void) {
    return current_memory_usage;
}

// Returns high-water mark
size_t inq_max_malloc_use(void) {
    return max_memory_usage;
}

// Reset high-water mark
void clear_max_malloc(void) {
    max_memory_usage = current_memory_usage;
}


void* tracked_realloc(void* ptr, size_t new_size) {
    if (!ptr) {
        // If ptr is NULL, behave like tracked_malloc
        return tracked_malloc(new_size);
    }

    if (new_size == 0) {
        // If new_size is 0, behave like tracked_free
        tracked_free(ptr);
        return NULL;
    }

    // Find the old size of the memory block
    size_t old_size = remove_allocation(ptr);

    // Reallocate the memory block
    void* new_ptr = realloc(ptr, new_size);
    if (!new_ptr) {
        // If realloc fails, reinsert the old allocation back into tracking
        add_allocation(ptr, old_size);
        fprintf(stderr, "Error: Memory reallocation failed.\n");
        return NULL;
    }

    // Update memory usage
    current_memory_usage += new_size - old_size;
    if (current_memory_usage > max_memory_usage) {
        max_memory_usage = current_memory_usage;
    }

    // Add the new allocation to the tracking system
    add_allocation(new_ptr, new_size);

    return new_ptr;
}


typedef enum {
    NC_UNSPECIFIED =  0,  /* ABSENT */
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
} NC_tag;

typedef struct {
    MPI_Offset  size;
    size_t      name_len; /* strlen(name), for faster string compare */
    char       *name;
} hdr_dim;


typedef struct hdr_dimarray {
    int            ndefined;      /* number of defined dimensions */
    // int            unlimited_id;  /* -1 for not defined, otherwise >= 0 */
    hdr_dim       **value;
} hdr_dimarray;

typedef struct {
    MPI_Offset nelems;   /* number of attribute elements */
    nc_type    xtype;    /* external NC data type of the attribute */
    size_t     name_len; /* strlen(name) for faster string compare */
    char      *name;     /* name of the attributes */
    void      *xvalue;   /* the actual data, in external representation */
} hdr_attr;


typedef struct hdr_attrarray {
    int            ndefined;  /* number of defined attributes */
    hdr_attr      **value;
} hdr_attrarray;


typedef struct {
    nc_type       xtype;   /* variable's external NC data type */
    size_t        name_len;/* strlen(name) for faster string compare */
    char         *name;    /* name of the variable */
    int           ndims;   /* number of dimensions */
    int          *dimids;  /* [ndims] array of dimension IDs */
    hdr_attrarray  attrs;   /* attribute array */
} hdr_var;


typedef struct hdr_vararray {
    int            ndefined;    /* number of defined variables */
    hdr_var       **value;
} hdr_vararray;


/* various file modes stored in flags */
struct hdr {
    MPI_Offset    xsz;      /* size occupied on the buffer */
    hdr_dimarray   dims;     /* dimensions defined */
    hdr_attrarray  attrs;    /* global attributes defined */
    hdr_vararray   vars;     /* variables defined */
};



typedef struct metabuffer {
    // MPI_Comm    comm;
    int         size;     /* allocated size of the buffer */
    char       *base;     /* beginning of read/write buffer */
    char       *pos;      /* current position in buffer */
    char       *end;      /* end position of buffer */
} metabuffer;


int
xlen_nc_type(nc_type xtype, int *size)
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
serialize_name(metabuffer *pbp,
                const char *name)
{
    size_t nchars = strlen(name);

    put_uint32((void**)(&pbp->pos), (uint32_t)nchars);

    return putn_text((void **)(&pbp->pos), (MPI_Offset)nchars, name);
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

    /* xlen_nc_type() returns the element size (unaligned) of
     * attrp->xtype attrp->xsz is the aligned total size of attribute values
     */
    xlen_nc_type(attrp->xtype, &xsz);
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



/*----< serialize_hdr() >----------------------------------------------*/
int
serialize_hdr(struct hdr *ncp, void *buf)
{
    int status;
    metabuffer putbuf;

    putbuf.pos           = buf;
    putbuf.base          = buf;
    putbuf.size          = ncp->xsz;

    /* copy dim_list */
    status = serialize_dimarray(&putbuf, &ncp->dims);
    if (status != NC_NOERR) return status;


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
    *name = (char *)tracked_malloc(nchars + 1);
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

        ncap->value = (hdr_dim **)tracked_malloc(ncap->ndefined * sizeof(hdr_dim *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }

        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_dim *)tracked_malloc(sizeof(hdr_dim));
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

    xlen_nc_type(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;

    attrp->xvalue = tracked_malloc(sz);
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
        ncap->value = (hdr_attr **)tracked_malloc(ncap->ndefined * sizeof(hdr_attr *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }
        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_attr *)tracked_malloc(sizeof(hdr_attr));
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
    varp->dimids = (int *)tracked_malloc(varp->ndims * sizeof(int));
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
        ncap->value = (hdr_var **)tracked_malloc(ncap->ndefined * sizeof(hdr_var *));
        if (ncap->value == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return -1;
        }

        for (int i = 0; i < ncap->ndefined; i++) {
            ncap->value[i] = (hdr_var *)tracked_malloc(sizeof(hdr_var));
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

int deserialize_hdr(struct hdr *ncp, void *buf, int buf_size) {

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
    assert((int)(getbuf.pos - getbuf.base) == getbuf.size);



    return 0;
}



void free_hdr_dim(hdr_dim *dim) {
    if (dim != NULL) {
        tracked_free(dim->name);
        tracked_free(dim);
    }
}

void free_hdr_dimarray(hdr_dimarray *dims) {
    if (dims != NULL) {
        for (int i = 0; i < dims->ndefined; i++) {
            
            free_hdr_dim(dims->value[i]);
        }
        tracked_free(dims->value);
        //tracked_free(dims);
    }
}

void free_hdr_attr(hdr_attr *attr) {
    if (attr != NULL) {
        tracked_free(attr->name);
        tracked_free(attr->xvalue);
    }
}

void free_hdr_attrarray(hdr_attrarray *attrs) {
    if (attrs != NULL) {
        if (attrs->value != NULL) {
            for (int i = 0; i < attrs->ndefined; i++) {
                free_hdr_attr(attrs->value[i]);
            }
            tracked_free(attrs->value);
            attrs->value = NULL;
            // tracked_free(attrs);
        }
    }
}

void free_hdr_var(hdr_var *var) {
    if (var != NULL) {
        tracked_free(var->name);
        tracked_free(var->dimids);

        free_hdr_attrarray(&(var->attrs));
        tracked_free(var);
    }
}

void free_hdr_vararray(hdr_vararray *vars) {
    if (vars != NULL) {
        for (int i = 0; i < vars->ndefined; i++) {
            free_hdr_var(vars->value[i]);
        }
        tracked_free(vars->value);
        // tracked_free(vars);
    }
}

void free_hdr(struct hdr *header) {
    if (header != NULL) {
        free_hdr_dimarray(&(header->dims));
        // free_hdr_attrarray(&(header->attrs));
        free_hdr_vararray(&(header->vars));
    }
}





/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("total maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   sum_size, (float)sum_size /1048576);
            printf("rank 0 maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            printf("rank 1 maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   malloc_size, (float)malloc_size /1048576);
        }
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}


static int
app_check_mem_usage(MPI_Comm comm)
{
    int err=0, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    malloc_size = inq_max_malloc_use();
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("total maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   sum_size, (float)sum_size /1048576);
            printf("rank 0 maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            printf("rank 1 maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   malloc_size, (float)malloc_size /1048576);
        }
    }
    return nerrs;
}

/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_crt_mem(MPI_Comm comm, int checkpoint)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("checkpoint %d: total current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   checkpoint, sum_size, (float)sum_size /1048576);
            printf("checkpoint %d: rank 0 current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   checkpoint, malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            printf("checkpoint %d: rank 1 current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   checkpoint, malloc_size, (float)malloc_size /1048576);
        }
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

/* check PnetCDF library internal memory usage */
static int
app_check_crt_mem(MPI_Comm comm, int checkpoint)
{
    int err=0, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    malloc_size = inq_malloc_use();
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("checkpoint %d: total current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   checkpoint, sum_size, (float)sum_size /1048576);
            printf("checkpoint %d: rank 0 current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   checkpoint, malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            printf("checkpoint %d: rank 1 current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   checkpoint, malloc_size, (float)malloc_size /1048576);
        }
    }
    return nerrs;
}





/* ---------------------------------- Generate Metadata ----------------------------------------*/

void generate_metadata(int rank, int nproc, struct hdr *file_info, int num_vars, int num_dims, int num_attrs, int dim_size, int att_size) {

    int ncid,  tot_num_dims, elem_sz, v_attrV_xsz, status;
    MPI_Offset start, count;
    file_info->dims.ndefined = 0;
    file_info->attrs.ndefined = 0;
    file_info->vars.ndefined = 0;
    file_info->xsz = 0;
    tot_num_dims = 0;

    file_info->xsz += 2 * sizeof(uint32_t); // NC_Variable and ndefined
    file_info->xsz += 2 * sizeof(uint32_t); // NC_Dimension and nelems
    // Calculate equal distribution of variables among processes
    num_vars = num_vars;
    int vars_per_process = num_vars / nproc;
    int remainder = num_vars % nproc;

    start = rank * vars_per_process + (rank < remainder ? rank : remainder);
    count = vars_per_process + (rank < remainder ? 1 : 0);

    file_info->vars.ndefined = count;
    file_info->vars.value = (hdr_var **)tracked_malloc(file_info->vars.ndefined * sizeof(hdr_var *));
    // Each process reads its subset of variables

    int dim_idx = 0;
    for (int i = start; i < start + count; ++i) {
        // Get variable information
        hdr_var *variable_info = (hdr_var *)tracked_malloc(sizeof(hdr_var));
        variable_info->ndims = num_dims;  // Initialize the number of dimensions
        variable_info->attrs.ndefined = 0;  // Initialize the number of attributes

        // Get variable information
        char var_name[NC_MAX_NAME + 1];
        sprintf(var_name, "process_%d_var_%d", rank, i - (int)start);
        variable_info->name_len = strlen(var_name);
        variable_info->name = (char *)tracked_malloc((variable_info->name_len + 1) * sizeof(char));
        strcpy(variable_info->name, var_name);

        variable_info->dimids = (int *)tracked_malloc(variable_info->ndims * sizeof(int));
        // Get dimension IDs
        for (int d = 0; d < variable_info->ndims; d++) {
            variable_info->dimids[d] = dim_idx;
            dim_idx++;
        }

        // Get varaible nc_type and attributes
        nc_type xtype = VAR_TYPE;
        
        variable_info->attrs.ndefined = num_attrs;
        variable_info->xtype = xtype;

        file_info->xsz += sizeof(uint32_t) + sizeof(char) * variable_info -> name_len; //var name
        file_info->xsz += sizeof(uint32_t); //xtype
        file_info->xsz += sizeof(uint32_t); //nelems of dim list
        file_info->xsz += sizeof(uint32_t) * variable_info ->ndims; // dimid list



        //Read and store dimension information 
        file_info->dims.ndefined += num_dims;
        if (tot_num_dims == 0) {
            file_info->dims.value = (hdr_dim **)tracked_malloc(file_info->dims.ndefined * sizeof(hdr_dim *));
        }else{
            file_info->dims.value = (hdr_dim **)tracked_realloc(file_info->dims.value, file_info->dims.ndefined * sizeof(hdr_dim *));
        } 

        for (int k = 0; k < num_dims; ++k) {
            hdr_dim *dimension_info = (hdr_dim *)tracked_malloc(sizeof(hdr_dim));
            int dimid = variable_info->dimids[k];
            // Get dimension name
            char dim_name[NC_MAX_NAME + 1];
            sprintf(dim_name, "process_%d_var_%d_dim_%d", rank, i, k);
            dimension_info->name_len = strlen(dim_name);
            dimension_info->name = (char *)tracked_malloc((dimension_info->name_len + 1) * sizeof(char));
            strcpy(dimension_info->name, dim_name);

            // Get dimension size

            dimension_info->size = dim_size;

            file_info->dims.value[k + tot_num_dims] = dimension_info;
            variable_info->dimids[k] = k + tot_num_dims; //overwriting previous global dim id to local dim id
            file_info->xsz += sizeof(uint32_t) + sizeof(char) * dimension_info -> name_len; // dim name
            file_info->xsz += sizeof(uint32_t); //size

        }
        tot_num_dims += num_dims;

        // Allocate memory for attributes
        file_info->xsz += 2 * sizeof(uint32_t); // NC_Attribute and ndefine
        variable_info->attrs.value = (hdr_attr **)tracked_malloc(num_attrs * sizeof(hdr_attr *));
        for (int j = 0; j < num_attrs; ++j) {
            variable_info->attrs.value[j] = (hdr_attr *)tracked_malloc(sizeof(hdr_attr));
            variable_info->attrs.value[j]->name_len = 0;  // Initialize attribute name length

            // Get attribute name 
            char att_name[NC_MAX_NAME + 1];
            sprintf(att_name, "process_%d_var_%d_attr_%d", rank, i, j);
            variable_info->attrs.value[j]->name_len = strlen(att_name);
            variable_info->attrs.value[j]->name = (char *)tracked_malloc(( strlen(att_name) + 1) * sizeof(char));
            strcpy(variable_info->attrs.value[j]->name, att_name);

            // Get attribute type and size
            nc_type attr_type = NC_INT;
            MPI_Offset attr_size = att_size;


            // Allocate memory for attribute value and read it
            variable_info->attrs.value[j]->xtype = attr_type;
            variable_info->attrs.value[j]->nelems = attr_size;
            xlen_nc_type(attr_type, &elem_sz);
            variable_info->attrs.value[j]->xvalue = tracked_malloc(attr_size * elem_sz);
            memset(variable_info->attrs.value[j]->xvalue, 0, sizeof(int) * attr_size);


            file_info->xsz += sizeof(uint32_t) + sizeof(char) * variable_info->attrs.value[j]->name_len; //attr name
            file_info->xsz += sizeof(uint32_t); // nc_type
            file_info->xsz += sizeof(uint32_t); // nelems
            status = xlen_nc_type(variable_info->attrs.value[j]->xtype, &v_attrV_xsz);
            file_info->xsz += variable_info->attrs.value[j]->nelems * v_attrV_xsz; // attr_value
        }

        // Add the variable information to the file_info structure
        file_info->vars.value[i - start] = variable_info;
    }

}
/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr(struct hdr *hdr_data, int ncid){
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int *dimid = (int *)tracked_malloc(ndims * sizeof(int));
    int i,j,k,nerrs=0;
    int err;

    for (i=0; i<ndims; i++){
        def_start_time = MPI_Wtime();
        err = ncmpi_def_dim(ncid, hdr_data->dims.value[i]->name,  hdr_data->dims.value[i]->size, &dimid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        // }
    }

    //define variables
    int nvars = hdr_data->vars.ndefined;
    int *varid = (int *)tracked_malloc(nvars * sizeof(int));
    int v_ndims, v_namelen, xtype, n_att;
    int *v_dimids;
    int att_namelen, att_xtype, att_nelems;

    for (i=0; i<nvars; i++){

        v_namelen =  hdr_data->vars.value[i]->name_len;
        xtype = hdr_data->vars.value[i]->xtype;

        v_ndims = hdr_data->vars.value[i]->ndims;
        v_dimids = (int *)tracked_malloc(v_ndims * sizeof(int));
        for(j=0; j<v_ndims; j++) v_dimids[j] = dimid[hdr_data->vars.value[i]->dimids[j]];
        def_start_time = MPI_Wtime();
        err = ncmpi_def_var(ncid, hdr_data->vars.value[i]->name, xtype, v_ndims,  v_dimids, &varid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        n_att = hdr_data->vars.value[i]->attrs.ndefined;

        for(k=0; k<n_att; k++){
            att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
            att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
            att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
            err = ncmpi_put_att(ncid, varid[i],  hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
            att_nelems, (const void *)hdr_data->vars.value[i]->attrs.value[k]->xvalue); ERR
        }
        tracked_free(v_dimids);

    }
    tracked_free(varid);
    tracked_free(dimid);
    return nerrs;
}

static int deserialize_all_hdr(struct hdr **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (struct hdr *)tracked_malloc(sizeof(struct hdr));
        deserialize_hdr(all_recv_hdr[i], all_collections_buffer + recv_displs[i], recvcounts[i]);
    }
    return 0;
}

int define_all_hdr(struct hdr **all_recv_hdr, int nproc, int ncid){
    for (int i=0; i< nproc; i++){
        struct hdr *hdr_data = all_recv_hdr[i];
        define_hdr(hdr_data, ncid);
    }
    return 0;
}

static int free_all_hdr(struct hdr **all_recv_hdr, int nproc){
    if (all_recv_hdr != NULL){
        for (int i=0; i< nproc; i++) free_hdr(all_recv_hdr[i]);
        tracked_free(all_recv_hdr);
    }
    return 0;
}




static void usage(const char *argv0)
{
    char *help =
    "Usage: %s [options]\n"
    "  [-q] <quite mode>\n"
    "  [-m] <memory tracking>\n"
    "  [-n] <num vars>      (default: %d)\n"
    "  [-d] <num dims per var>      (default: %d)\n"
    "  [-a] <num attrs per var>     (default: %d)\n"
    "  [-l] <dim size>      (default: %d)\n"
    "  [-s] <attr size>     (default: %d)\n"
    "  [-t] <hash table size>   (default: %d)\n"
    "  [-h]                 Show this help message\n"
    "  [filename]: output netCDF file name (default %s)\n\n";

    fprintf(stderr,
            help,
            argv0,
            NUM_VARS,
            NUM_DIMS,
            NUM_ATTRS,
            DIM_SIZE,
            ATTR_SIZE,
            HASH_SIZE,
            OUTPUT_NAME);

}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, i, nerrs=0;
    char filename[256];
    double end_to_end_time, mpi_time, io_time, enddef_time, close_time, max_time, min_time;
    double start_time, start_time1, end_time1, end_time2, end_time3, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    struct hdr local_hdr;


    int verbose=1, mem_track=0, num_vars=NUM_VARS, num_dims_per_var=NUM_DIMS, num_attrs_per_var=NUM_ATTRS, dim_size=DIM_SIZE, attr_size=ATTR_SIZE, hash_size=HASH_SIZE;

    while ((i = getopt(argc, argv, "hqma:d:l:s:t:n:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'm': mem_track = 1;
                      break;
            case 'n': num_vars = atoi(optarg);
                      break;
            case 'a': num_attrs_per_var = atoi(optarg);
                      break;
            case 'd': num_dims_per_var = atoi(optarg);
                      break;
            case 'l': dim_size = atoi(optarg);
                      break;
            case 's': attr_size = atoi(optarg);
                      break;
            case 't': hash_size = atoi(optarg);
                      break;
            case 'h': 
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, OUTPUT_NAME);
    else                      snprintf(filename, 256, "%s", argv[optind]);

    generate_metadata(rank, nproc, &local_hdr, num_vars, num_dims_per_var, num_attrs_per_var, dim_size, attr_size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = start_time1 = MPI_Wtime();
    char* send_buffer = (char*) tracked_malloc(local_hdr.xsz);
    status = serialize_hdr(&local_hdr, send_buffer);

    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) tracked_malloc(nproc * sizeof(MPI_Offset));
    MPI_Allgather(&local_hdr.xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

    // Calculate displacements for the second phase
    int* recv_displs = (int*) tracked_malloc(nproc * sizeof(int));
    int total_recv_size, min_size, max_size;
    total_recv_size = min_size = max_size = all_collection_sizes[0];
    recv_displs[0] = 0;


    for (int i = 1; i < nproc; ++i) {
        recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
        total_recv_size += all_collection_sizes[i];
        if(all_collection_sizes[i] > max_size){
            max_size = all_collection_sizes[i];
        }
        if(all_collection_sizes[i] < min_size){
            min_size = all_collection_sizes[i];
        }
    }
    double total_recv_size_MB = total_recv_size / (1024.0 * 1024.0);
    double min_size_MB = min_size / (1024.0 * 1024.0);
    double max_size_MB = max_size / (1024.0 * 1024.0);

    char* all_collections_buffer = (char*) tracked_malloc(total_recv_size);
    int* recvcounts =  (int*)tracked_malloc(nproc * sizeof(int));
    for (int i = 0; i < nproc; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    // Phase 2: Communicate the actual header data
    // Before MPI_Allgatherv
    MPI_Allgatherv(send_buffer, local_hdr.xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);
    // Deserialize the received data and print if rank is 0
    
    int ncid, cmode;

    cmode = NC_64BIT_DATA | NC_CLOBBER;

    struct hdr **all_recv_hdr = (struct hdr **)tracked_malloc(nproc * sizeof(struct hdr*));
    deserialize_all_hdr(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    char hash_size_str[64];
    snprintf(hash_size_str, sizeof(hash_size_str), "%d", hash_size);
    MPI_Info_set(info, "nc_hash_size_dim", hash_size_str);
    MPI_Info_set(info, "nc_hash_size_var", hash_size_str);
    if (rank == 0 && verbose){ 
        printf("Hash table size for dim: %d\n", hash_size);
        printf("Hash table size for var: %d\n", hash_size);
    }
    if (mem_track){
        app_check_crt_mem(MPI_COMM_WORLD, 0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid); ERR
    MPI_Info_free(&info);

    define_all_hdr(all_recv_hdr, nproc, ncid);
    if (mem_track){
        app_check_crt_mem(MPI_COMM_WORLD, 1);
        pnetcdf_check_crt_mem(MPI_COMM_WORLD, 1);
    }




    io_time = MPI_Wtime() - end_time1;

    free_all_hdr(all_recv_hdr, nproc);
    if (mem_track){
    app_check_crt_mem(MPI_COMM_WORLD, 2);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 2);
    }
    end_time2 = MPI_Wtime();

    err = ncmpi_enddef(ncid); ERR

    end_time3 = MPI_Wtime();
    if (mem_track){
    app_check_crt_mem(MPI_COMM_WORLD, 3);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 3);
    }
    enddef_time = end_time3 - end_time2;


    // Clean up
    free_hdr(&local_hdr);
    tracked_free(send_buffer);
    tracked_free(all_collections_buffer);
    tracked_free(all_collection_sizes);
    tracked_free(recv_displs);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time3 = MPI_Wtime();
    if (mem_track){
    app_check_crt_mem(MPI_COMM_WORLD, 4);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 4);
    }
    err = ncmpi_close(ncid); ERR
    end_time = MPI_Wtime();
    close_time = end_time - end_time3;
    end_to_end_time = end_time - start_time;
    mpi_time = end_time1 - start_time1;


    double times[6] = {end_to_end_time, mpi_time, io_time, enddef_time, total_def_time, close_time};
    char *names[6] = {"end-end", "mpi-phase", "write", "enddef", "def_dim/var", "close"};
    double max_times[6], min_times[6];


    MPI_Reduce(&times[0], &max_times[0], 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 6, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 6; i++) {
        if (rank == 0 && verbose) {
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
        }
    }
    if (mem_track){
    app_check_crt_mem(MPI_COMM_WORLD, 5);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 5);
    pnetcdf_check_mem_usage(MPI_COMM_WORLD);
    app_check_mem_usage(MPI_COMM_WORLD);
    }
    free_allocation_struct();
    MPI_Finalize();
    return 0;
}