/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program shows how to use the one-file-per-thread I/O operation. In this
 * example, each MPI process creates 6 POSIX threads and each thread does the
 * followings (one unique file per thread):
 *   1. creates a unique new file,
 *   2. writes 2 records to a record variable
 *   3. writes a fixed-size variable,
 *   4. closes the file,
 *   5. re-open the file,
 *   6. reads the record variable and checks contents,
 *   7. reads the fixed-size variable and checks contents,
 *   8. closes file.
 *
 * To compile:
 *    % mpicc -O2 pthread.c -o pthread -lpnetcdf -lpthread
 *
 * Example commands for MPI run on 4 MPI processes.
 *    % mpiexec -n 4 ./pthread testfile.nc
 *
 * This example run will create 24 files.
 *    % ls -lgG .
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.0
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.1
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.4
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.2
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.3
 *    -rw------- 1    1072 Jul 21 14:33 testfile.nc.5
 *    ...
 *
 * The file header of all output files looks the same. One example from command
 * ncmpidump is given below.
 *    % ncdump testfile.nc.2
 *    netcdf testfile.nc {
 *    dimensions:
 *            time = UNLIMITED ; // (3 currently)
 *            X = 4 ;
 *            Y = 5 ;
 *    variables:
 *            int ivar(time, X) ;
 *            double dvar(Y, X) ;
 *    data:
 *
 *     ivar =
 *      2, 2, 2, 2,
 *      0, 0, 0, 0,
 *      2, 2, 2, 2 ;
 *
 *     dvar =
 *      2, 2, 2, 2,
 *      2, 2, 2, 2,
 *      2, 2, 2, 2,
 *      2, 2, 2, 2,
 *      2, 2, 2, 2 ;
 *    }
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* _POSIX_BARRIERS, getopt() */
#include <pthread.h>

#include <mpi.h>
#include <pnetcdf.h>

#define NTHREADS 6
#define NY 5
#define NX 4

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

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

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name base\n";
    fprintf(stderr, help, argv0);
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
        if (rank == 0 && verbose)
            printf("maximum heap memory allocated by PnetCDF internally is %lld bytes\n",
                   sum_size);

        /* check if there is any PnetCDF internal malloc residue */
        err = ncmpi_inq_malloc_size(&malloc_size);
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

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
    int i, id, nprocs, err, nerrs=0, ncid, *ret, dimid[2], varid[2];
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
    err = ncmpi_create(MPI_COMM_SELF, filename, NC_CLOBBER, info, &ncid); ERR
    MPI_Info_free(&info);

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[1]); ERR
    /* define a record variable ivar of integer type */
    err = ncmpi_def_var(ncid, "ivar", NC_INT, 2, dimid, &varid[0]); ERR
    /* define a fixed-size variable dvar of double type */
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "dvar", NC_DOUBLE, 2, dimid, &varid[1]); ERR
    err = ncmpi_enddef(ncid); ERR
    /* now we are in data mode */

    /* write a record to the record variable */
    start[0] = 0; /* first record */
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, ibuf); ERR

    /* write another record to the record variable */
    start[0] = 2; /* third record */
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_int_all(ncid, varid[0], start, count, ibuf); ERR

    /* write to the fixed-size variable */
    err = ncmpi_put_var_double_all(ncid, varid[1], dbuf); ERR

    err = ncmpi_close(ncid); ERR

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

    err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); ERR
    err = ncmpi_inq_varid(ncid, "ivar", &varid[0]); ERR
    err = ncmpi_inq_varid(ncid, "dvar", &varid[1]); ERR

    /* read the first record of the record variable */
    for (i=0; i<NX; i++) ibuf[i] = -1;
    start[0] = 0;
    start[1] = 0;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, ibuf); ERR
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
    err = ncmpi_get_vara_int_all(ncid, varid[0], start, count, ibuf); ERR
    for (i=0; i<NX; i++) {
        if (ibuf[i] != id) {
            printf("Error at %s line %d: expect ibuf[%d]=%d but got %d\n",
            __FILE__, __LINE__, i, id, ibuf[i]);
            nerrs++;
            break;
        }
    }

    /* read the fixed-size variable */
    err = ncmpi_get_var_double_all(ncid, varid[1], dbuf); ERR
    for (i=0; i<NY*NX; i++) {
        if (dbuf[i] != (double)id) {
            printf("Error at %s line %d: expect ibuf[%d]=%d but got %f\n",
            __FILE__, __LINE__, i, id, dbuf[i]);
            nerrs++;
            break;
        }
    }
    err = ncmpi_close(ncid); ERR

    free(ibuf);
    free(dbuf);

    /* return number of errors encountered */
    ret = (int*)malloc(sizeof(int));
    *ret = nerrs;

    return ret; /* same as pthread_exit(ret); */
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
    extern int optind;
    char filename[256];
    int  i, nerrs=0, rank, providedT;
    pthread_t threads[NTHREADS];

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedT);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;

    if(rank == 0 && verbose) {
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
    if (providedT != MPI_THREAD_MULTIPLE) {
        if (!rank)
            printf("\nWarning: MPI provided thread support level is less than MPI_THREAD_MULTIPLE ---- skip this test\n");
        MPI_Finalize();
        return 0;
    }

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

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

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

