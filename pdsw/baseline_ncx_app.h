#include <stddef.h> 
#include <mpi.h>
#include <pnetcdf.h>


#define NC_NOERR        0  


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


// Function prototypes
int xlen_nc_type_meta(nc_type xtype, int *size);
static int putn_text(void **xpp, MPI_Offset nelems, const char *tp);
static int put_uint32(void **xpp, unsigned int ip);
static int serialize_dim(metabuffer *pbp, const hdr_dim *dimp);
static int serialize_name(metabuffer *pbp, const char *name);
static int serialize_dimarray(metabuffer *pbp, const hdr_dimarray *ncap);
static int serialize_attrV(metabuffer *pbp, const hdr_attr *attrp);
static int serialize_attr(metabuffer *pbp, const hdr_attr *attrp);
static int serialize_attrarray(metabuffer *pbp, const hdr_attrarray *ncap);
static int serialize_var(metabuffer *pbp, const hdr_var *varp);
static int serialize_vararray(metabuffer *pbp, const hdr_vararray *ncap);
int serialize_hdr_meta(struct hdr *ncp, void *buf);


// Deserialization functions
int deserialize_hdr_meta(struct hdr *ncp, void *buf, int buf_size);
static int getn_text(void **xpp, MPI_Offset nelems, char *tp);
static int get_uint32(void **xpp, unsigned int *ip);
static int deserialize_nc_type(metabuffer *gbp, nc_type *xtypep);
static int deserialize_name(metabuffer *gbp, char **name);
static int deserialize_dim(metabuffer *gbp, hdr_dim *dimp);
static int deserialize_dimarray(metabuffer *gbp, hdr_dimarray *ncap);
static int deserialize_attrV(metabuffer *gbp, hdr_attr *attrp);
static int deserialize_attr(metabuffer *gbp, hdr_attr *attrp);
static int deserialize_attrarray(metabuffer *gbp, hdr_attrarray *ncap);
static int deserialize_var(metabuffer *gbp, hdr_var *varp);
static int deserialize_vararray(metabuffer *gbp, hdr_vararray *ncap);

// Read and distribute metadata

int read_metadata_from_file(const char* filename, struct hdr *recv_hdr);
int distribute_metadata(int rank, int nproc, struct hdr *all_hdr, struct hdr *local_hdr);

/* Function prototypes */
void free_hdr_dim_meta(hdr_dim *dim);
void free_hdr_dimarray_meta(hdr_dimarray *dims);
void free_hdr_attr_meta(hdr_attr *attr);
void free_hdr_attrarray_meta(hdr_attrarray *attrs);
void free_hdr_var_meta(hdr_var *var);
void free_hdr_vararray_meta(hdr_vararray *vars);
void free_hdr_meta(struct hdr *header);
