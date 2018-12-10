#define NC_ZIP_DRIVER_DUMMY 0

struct NCZIP_driver {
    int (*init)(MPI_Info);
    int (*finalize)(MPI_Comm, const char*, int, int, MPI_Info, void**);
    int (*compress)(void*, MP_Offset, void*, MP_Offset*, int, int*, MPI_Datatype);
    int (*decompress)(void*, MP_Offset, void*, MP_Offset*, int, int*, MPI_Datatype);
};

typedef struct NCZIP_driver NCZIP_driver;

extern NCZIP_driver* nczip_dummy_inq_driver(void);