#ifndef NCCHK_FILTER_DRIVER_H
#define NCCHK_FILTER_DRIVER_H

#include <mpi.h>

struct NCCHK_filter {
    int (*init)(MPI_Info);
    int (*finalize)();
    int (*inq_cpsize)(void*, int, int*, int, int*, MPI_Datatype);
    int (*compress)(void*, int, void*, int*, int, int*, MPI_Datatype);
    int (*compress_alloc)(void*, int, void**, int*, int, int*, MPI_Datatype);
    int (*inq_dcsize)(void*, int, int*, int, int*, MPI_Datatype);
    int (*decompress)(void*, int, void*, int*, int, int*, MPI_Datatype);
    int (*decompress_alloc)(void*, int, void**, int*, int, int*, MPI_Datatype);
};

typedef struct NCCHK_filter NCCHK_filter;

extern NCCHK_filter* ncchk_dummy_inq_driver(void);

#if ENABLE_ZLIB
extern NCCHK_filter* ncchk_zlib_inq_driver(void);
#endif

#if ENABLE_SZ
extern NCCHK_filter* ncchk_sz_inq_driver(void);
#endif

#endif