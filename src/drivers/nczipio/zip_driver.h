#define NC_ZIP_DRIVER_DUMMY 0


struct NCZIP_driver {
    int (*init)(MPI_Info);
    int (*finalize)();
    int (*compress)(void*, int, void*, int*, int, int*, MPI_Datatype);
    int (*decompress)(void*, int, void*, int*, int, int*, MPI_Datatype);
};

typedef struct NCZIP_driver NCZIP_driver;

extern NCZIP_driver* nczip_dummy_inq_driver(void);