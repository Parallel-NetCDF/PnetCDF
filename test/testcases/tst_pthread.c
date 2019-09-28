/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program tests thread-safe capability. Each MPI process creates 6
 * threads and each thread does the followings (one unique file per thread):
 * 1. creates a unique new file,
 * 2. writes 2 records to a record variable
 * 3. writes a fixed-size variable,
 * 4. closes the file,
 * 5. open a different file created by another thread,
 * 6. reads the record variable and check contents,
 * 7. reads the fixed-size variable and check contents,
 * 8. closes the file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <unistd.h> /* _POSIX_BARRIERS, unlink() */

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#ifdef ENABLE_THREAD_SAFE
#include <pthread.h>

#define NTHREADS 6
#define NY 5
#define NX 4

#if !defined(_POSIX_BARRIERS) || _POSIX_BARRIERS <= 0
/* According to opengroup.org, barriers are defined in the optional part of
 * POSIX standard. For example, Mac OSX does not have pthread_barrier. If
 * barriers were implemented, the _POSIX_BARRIERS macro is defined as a
 * positive number.
 */

#include <errno.h>

typedef int pthread_barrierattr_t;
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int numThreads;
} pthread_barrier_t;

static int pthread_barrier_init(pthread_barrier_t           *barrier,
                                const pthread_barrierattr_t *attr,
                                unsigned int                 count)
{
    if (count == 0) {
        errno = EINVAL;
        return -1;
    }

    if (pthread_mutex_init(&barrier->mutex, 0) < 0)
        return -1;

    if (pthread_cond_init(&barrier->cond, 0) < 0) {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }
    pthread_mutex_lock(&barrier->mutex);
    barrier->numThreads = count;
    barrier->count = 0;
    pthread_mutex_unlock(&barrier->mutex);

    return 0;
}

static int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

static int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    int ret;
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if (barrier->count >= barrier->numThreads) {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        ret = 1;
    } else {
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
        ret = 0;
    }
    pthread_mutex_unlock(&barrier->mutex);
    return ret;
}
#endif

/* pthread barrier object */
static pthread_barrier_t barr;

typedef struct {
    int  id;         /* globally unique thread ID */
    char fname[256]; /* output file name base */
} thread_arg;

/*----< thread_func() >------------------------------------------------------*/
static
void* thread_func(void *arg)
{
    char filename[256];
    int i, id, nprocs, cmode, err, nerrs=0, ncid, *ret, dimid[2], varid[2];
    int *ibuf;
    double *dbuf;
    MPI_Offset start[2], count[2];
    MPI_Info info;

    /* make a unique file name for each thread */
    id = ((thread_arg*)arg)->id;
    sprintf(filename, "%s.%d", ((thread_arg*)arg)->fname, id);

    /* allocate I/O buffers and initialize their contents */
    ibuf = (int*)    malloc(NY * NX * sizeof(int));
    dbuf = (double*) malloc(NY * NX * sizeof(double));
    for (i=0; i<NY*NX; i++) {
        ibuf[i] = id;
        dbuf[i] = 1.0 * id;
    }

    /* set an MPI-IO hint to disable file offset alignment for fixed-size
     * variables */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_var_align_size", "1");

    /* create a file, clobber it if already exits */
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_SELF, filename, cmode, info, &ncid); CHECK_ERR
    MPI_Info_free(&info);

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); CHECK_ERR
    /* define a record variable ivar of integer type */
    err = ncmpi_def_var(ncid, "ivar", NC_INT, 2, dimid, &varid[0]); CHECK_ERR
    /* define a fixed-size variable dvar of double type */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "dvar", NC_DOUBLE, 2, dimid, &varid[1]); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR
    /* now we are in data mode */

    /* write a record to the record variable */
    start[0] = 0; /* first record */
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, ibuf); CHECK_ERR

    /* write another record to the record variable */
    start[0] = 2; /* third record */
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, ibuf); CHECK_ERR

    /* write to the fixed-size variable */
    err = ncmpi_put_var_double_all(ncid, varid[1], dbuf); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    /* synchronize all processes (only one thread per process participates) */
    if (id % NTHREADS == 0) MPI_Barrier(MPI_COMM_WORLD);

    /* synchronize all threads within each process to ensure all threads to
     * finish their file writes */
    pthread_barrier_wait(&barr);

    /* each thread opens a different file (round-robin shift), reads variables
     * and check contents */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    id = (id + 1) % (nprocs * NTHREADS);
    sprintf(filename, "%s.%d", ((thread_arg*)arg)->fname, id);

    err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "ivar", &varid[0]); CHECK_ERR
    err = ncmpi_inq_varid(ncid, "dvar", &varid[1]); CHECK_ERR

    /* read the first record of the record variable */
    for (i=0; i<NX; i++) ibuf[i] = -1;
    start[0] = 0;
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, ibuf); CHECK_ERR
    for (i=0; i<NX; i++) {
        if (ibuf[i] != id) {
            printf("Error at %s line %d: expect ibuf[%d]=%d but got %d\n",
            __FILE__, __LINE__, i, id, ibuf[i]);
            nerrs++;
            break;
        }
    }

    /* read the 3rd record of the record variable */
    for (i=0; i<NX; i++) ibuf[i] = -1;
    start[0] = 2;
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, ibuf); CHECK_ERR
    for (i=0; i<NX; i++) {
        if (ibuf[i] != id) {
            printf("Error at %s line %d: expect ibuf[%d]=%d but got %d\n",
            __FILE__, __LINE__, i, id, ibuf[i]);
            nerrs++;
            break;
        }
    }

    /* read the fixed-size variable */
    err = ncmpi_get_var_double_all(ncid, varid[1], dbuf); CHECK_ERR
    for (i=0; i<NY*NX; i++) {
        if (dbuf[i] != (double)id) {
            printf("Error at %s line %d: expect ibuf[%d]=%d but got %f\n",
            __FILE__, __LINE__, i, id, dbuf[i]);
            nerrs++;
            break;
        }
    }
    err = ncmpi_close(ncid); CHECK_ERR

    free(ibuf);
    free(dbuf);

    // unlink(filename);

    /* return number of errors encountered */
    ret = (int*)malloc(sizeof(int));
    *ret = nerrs;
    return ret; /* same as pthread_exit(ret); */
}
#endif

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
    char filename[256];
    int  i, err, nerrs=0, rank, providedT;

#ifdef ENABLE_THREAD_SAFE
    pthread_t threads[NTHREADS];
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedT);
#else
    MPI_Init(&argc, &argv);
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for thread safety ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

#ifdef ENABLE_THREAD_SAFE
#ifdef DEBUG
    if (rank == 0) {
        switch (providedT) {
            case MPI_THREAD_SINGLE:      printf("Support MPI_THREAD_SINGLE\n");
                                         break;
            case MPI_THREAD_FUNNELED:    printf("Support MPI_THREAD_FUNNELED\n");
                                         break;
            case MPI_THREAD_SERIALIZED:  printf("Support MPI_THREAD_SERIALIZED\n");
                                         break;
            case MPI_THREAD_MULTIPLE:    printf("Support MPI_THREAD_MULTIPLE\n");
                                         break;
            default: printf("Error MPI_Init_thread()\n"); break;
        }
    }
#endif
    if (providedT != MPI_THREAD_MULTIPLE) {
        if (!rank) {
            char fname[512];
            printf("\nWarning: MPI provided thread support level is less than MPI_THREAD_MULTIPLE ---- skip this test\n");
            for (i=0; i<NTHREADS; i++) { /* create dummy files for ncvalidator to check */
                int ncid;
                sprintf(fname, "%s.%d", filename, i);
                err = ncmpi_create(MPI_COMM_SELF, fname, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR
                err = ncmpi_close(ncid); CHECK_ERR
            }
        }
        MPI_Finalize();
        return 0;
    }

    /* initialize thread barrier */
    pthread_barrier_init(&barr, NULL, NTHREADS);

    /* create threads, each calls thread_func() */
    for (i=0; i<NTHREADS; i++) {
        thread_arg t_arg[NTHREADS]; /* must be unique to each thread */
        t_arg[i].id = i + rank * NTHREADS;
        sprintf(t_arg[i].fname, "%s",filename);
        if (pthread_create(&threads[i], NULL, thread_func, &t_arg[i])) {
            fprintf(stderr, "Error in %s line %d creating thread %d\n",
                    __FILE__, __LINE__, i);
            nerrs++;
        }
    }

    /* wait for all threads to finish */
    for (i=0; i<NTHREADS; i++) {
        void *ret;
        if (pthread_join(threads[i], (void**)&ret)) {
            fprintf(stderr, "Error in %s line %d joining thread %d\n",
                    __FILE__, __LINE__, i);
        }
        nerrs += *(int*)ret;
        free(ret);
    }

    pthread_barrier_destroy(&barr);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }
#else
    if (rank == 0) printf(SKIP_STR);
#endif

    MPI_Finalize();
    return (nerrs > 0);
}
