/******************************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>

#include <mpi.h>

#define NTHREADS 6

#define ERRNO_HANDLE(errno) {                                   \
    if (errno != 0) {                                           \
        fprintf(stderr,"Error: %s (file=%s line=%d func=%s)\n", \
                strerror(errno),__FILE__,__LINE__,__func__);    \
        assert(err == 0);                                       \
    }                                                           \
}

#define ERR_HANDLER(err, err_msg) { \
    if (err != MPI_SUCCESS) { \
        int rank, errorclass, errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
        MPI_Error_string(err, errorString, &errorStringLen); \
        printf("rank %d: MPI error (%s) : %s\n", rank, err_msg, errorString); \
        assert(err == MPI_SUCCESS); \
    } \
}

typedef struct {
    int  id;         /* globally unique thread ID */
    char fname[256]; /* output file name base */
} thread_arg;

static
void* thread_func(void *arg)
{
    char filename[512];
    int id, err;
    MPI_File fh;
    MPI_Info info;

    /* make a unique file name for each thread */
    id = ((thread_arg*)arg)->id;
    sprintf(filename, "%s.%d", ((thread_arg*)arg)->fname, id);

    err = MPI_Info_create(&info);
    ERR_HANDLER(err, "MPI_Info_create")
    err = MPI_Info_set(info, "romio_ds_write", "disable");
    ERR_HANDLER(err, "MPI_Info_set")

    err = MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
                        info, &fh);
    ERR_HANDLER(err, "MPI_File_open")

    err = MPI_Info_free(&info);
    ERR_HANDLER(err, "MPI_Info_free")

    err = MPI_File_close(&fh);
    ERR_HANDLER(err, "MPI_File_close")

    return NULL;
}

static void thread_loop(int n)
{
    char *filename = "ufs:testfile";
    int i, err, rank;
    pthread_t threads[NTHREADS];

    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ERR_HANDLER(err, "MPI_Comm_rank")

    // err = MPI_Barrier(MPI_COMM_WORLD);
    // ERR_HANDLER(err, "MPI_Barrier")

    for (i=0; i<NTHREADS; i++) {
        thread_arg t_arg[NTHREADS]; /* must be unique to each thread */
        t_arg[i].id = n + i + rank * NTHREADS;
        sprintf(t_arg[i].fname, "%s",filename);
        err = pthread_create(&threads[i], NULL, thread_func, &t_arg[i]);
        ERRNO_HANDLE(err)
    }

    /* wait for all threads to finish */
    for (i=0; i<NTHREADS; i++) {
        void *ret;
        err = pthread_join(threads[i], (void**)&ret);
        ERRNO_HANDLE(err)
    }
}

int main(int argc, char *argv[])
{
    int  i, err, rank, nprocs, providedT;

    err = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedT);
    ERR_HANDLER(err, "MPI_Init_thread")

    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ERR_HANDLER(err, "MPI_Comm_rank")
    err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ERR_HANDLER(err, "MPI_Comm_size")

    if (providedT != MPI_THREAD_MULTIPLE) {
        if (!rank)
            printf("\nWarning: MPI provided thread support level is less than MPI_THREAD_MULTIPLE ---- skip this test\n");
        MPI_Finalize();
        exit(1);
    }

    for (i=0; i<6; i++)
        thread_loop(i*nprocs*NTHREADS);

    MPI_Finalize();
    exit(0);
}

