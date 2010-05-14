/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <string.h>
#include <unistd.h>

/* need this for SIZEOF_INT so we get the correct printf output */
#include "ncconfig.h"

#ifdef SIZEOF_INT
# if SIZEOF_INT == 4
#  define lld(x) (x)
# elif  SIZEOF_INT == 8
#  define lld(x) (long long)(x)
# endif
#endif

/* The file name is taken as a command-line argument. */

/* Measures the I/O bandwidth for writing/reading a 3D
   block-distributed array to a file corresponding to the global array
   in row-major (C) order.
   Note that the file access pattern is noncontiguous.

   Array size 128^3. For other array sizes, change array_of_gsizes below.*/

void* xmalloc(size_t size) {
    void *buf;
    buf = malloc(size);
    if (buf == NULL) {
        printf("Error: malloc of size=%d failed\n",size);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
    return buf;
}

#define TEST_HANDLE_ERR(status)                                   \
{                                                                 \
    if ((status) != NC_NOERR)                                     \
        printf("%s at line %d of %s\n", ncmpi_strerror((status)), \
                __LINE__, __FILE__);                              \
}


int main(int argc, char **argv)
{
    int i, j, array_of_gsizes[3],array_of_distribs[3];
    int order, nprocs, len, **buf, mynod;
    MPI_Offset bufcount;
    int array_of_dargs[3], array_of_psizes[3];
    int status;
    MPI_Offset sizes[3], array_of_starts[3];
    double write_time, *new_write_tim;
    double start_time, open_time, def_time, run_time;
    double *new_open_tim, *new_def_tim, *new_run_tim;
    char *pathname, filename[50];
    char dimname[20], varname[20];
    int ncid, dimids[3], rank_dim[3], *varid;
    MPI_Info info;
    MPI_Offset **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    int ndims = 3;
    int nvars = 10;
    int k, k_loop;
    MPI_Datatype *datatype_list;
    int length;
    int mvar_flag = 0;
    int *array_of_statuses, *array_of_requests;
    int unlimit_flag;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynod);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc < 6) {
        if (!mynod)
            printf("Usage: <array length> <number variables> <number loops> <API type> <unlimit type> <-fname pathname>\n");
        MPI_Finalize();
        return 1;
    }

    length       = atoi(argv[1]);
    nvars        = atoi(argv[2]);
    k_loop       = atoi(argv[3]);
    mvar_flag    = atoi(argv[4]);
    unlimit_flag = atoi(argv[5]);

    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

/* process 0 takes the file name as a command-line argument and
   broadcasts it to other processes */
    if (!mynod) {
        i = 1;
        while ((i < argc) && strcmp("-fname", *argv)) {
            i++;
            argv++;
        }
        if (i >= argc) {
            fprintf(stderr, "\n*#  Usage: coll_perf -fname pathname\n\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        argv++;
        len = strlen(*argv);
        pathname = (char *) xmalloc(len+1);

        strcpy(pathname, *argv);
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(pathname, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        pathname = (char *) xmalloc(len+1);
        MPI_Bcast(pathname, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    order = MPI_ORDER_C;

    buf = (int **)xmalloc(nvars*sizeof(int*));
    varid = (int *)xmalloc(nvars*sizeof(int));
    bufcount_list = (MPI_Offset *)xmalloc(nvars*sizeof(MPI_Offset));
    starts_list = (MPI_Offset **)xmalloc(nvars*sizeof(MPI_Offset *));
    count_list = (MPI_Offset **)xmalloc(nvars*sizeof(MPI_Offset *));
    datatype_list = (MPI_Datatype*)xmalloc(nvars*sizeof(MPI_Datatype));

    array_of_requests = (int *)xmalloc(2*nvars*sizeof(int));
    array_of_statuses = array_of_requests + nvars;

    new_open_tim  = (double *)xmalloc(4*k_loop*sizeof(double));
    new_def_tim   = new_open_tim  + k_loop;
    new_write_tim = new_def_tim   + k_loop;
    new_run_tim   = new_write_tim + k_loop;

    for (i=0; i<nvars; i++) {
        starts_list[i] = (MPI_Offset *)xmalloc(ndims*sizeof(MPI_Offset));
        count_list[i] = (MPI_Offset *)xmalloc(ndims*sizeof(MPI_Offset));
    }

    array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;

    array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;

    bufcount = 1;
    for (i=0; i<ndims; i++) {
        array_of_psizes[i] = 0;
        sizes[i] = length;
        bufcount *= length;
    }
    MPI_Dims_create(nprocs, ndims, array_of_psizes);

   /* subarray in each process is len x len x len */
    for (i=0; i<ndims; i++)
        array_of_gsizes[i] = length * array_of_psizes[i];

    /* mynd's process rank in each dimension (in MPI_ORDER_C) */
    rank_dim[2] =  mynod %  array_of_psizes[2];
    rank_dim[1] = (mynod /  array_of_psizes[2]) % array_of_psizes[1];
    rank_dim[0] =  mynod / (array_of_psizes[2]  * array_of_psizes[1]);

    /* starting coordinates of the subarray in each dimension */
    for (i=0; i<ndims; i++)
        array_of_starts[i] = length * rank_dim[i];

    /* mput */
    for (i=0; i<nvars; i++) {
        for (j=0; j<ndims; j++) {
           starts_list[i][j] = array_of_starts[j];
           count_list[i][j]  = length;
        }
        bufcount_list[i] = bufcount;
        datatype_list[i] = MPI_INT;
    }

    srand(mynod+1);
    for (i=0; i<nvars;i++){
        buf[i] = (int *) xmalloc(bufcount * sizeof(int));
        for (j=0; j<bufcount; j++)
            buf[i][j]=rand();
//          buf[i][j]=mynod+1;
//          buf[i][j]= mynod + 1 + 32768*i;
    }

    MPI_Info_create(&info);
/*
    MPI_Info_set(info, "cb_buffer_size", "1024");
    MPI_Info_set(info, "cb_buffer_size", "16777216");
    MPI_Info_set(info, "romio_no_indep_rw", "true");
    MPI_Info_set(info, "romio_cb_write", "true");
*/

    for (k=0; k<k_loop; k++){
        sprintf(filename, "%s.%d.%d.%d.%d.nc", pathname, length, nvars,
                mvar_flag, k);
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
        status = ncmpi_create(MPI_COMM_WORLD, filename,
                              NC_CLOBBER|NC_64BIT_OFFSET, info, &ncid);
        TEST_HANDLE_ERR(status);
        /* define dimensions */
        if (unlimit_flag == 1) {
            sprintf(dimname, "dim_%d", 0);
            status = ncmpi_def_dim(ncid, dimname, NC_UNLIMITED, &dimids[0]);
            TEST_HANDLE_ERR(status);
            for (i=1; i<ndims; i++){
                sprintf(dimname, "dim_%d", i);
                status = ncmpi_def_dim(ncid, dimname, array_of_gsizes[i],
                                       &dimids[i]);
                TEST_HANDLE_ERR(status);
            }
        } else {
            for (i=0; i<ndims; i++){
                sprintf(dimname, "dim_%d", i);
                status = ncmpi_def_dim(ncid, dimname, array_of_gsizes[i],
                                       &dimids[i]);
                TEST_HANDLE_ERR(status);
            }
        }

        /* define variables */
        for (i=0; i<nvars; i++){
            sprintf(varname, "var_%d", i);
            status = ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids,
                                   &varid[i]);
            TEST_HANDLE_ERR(status);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        open_time = MPI_Wtime()-start_time;

        status = ncmpi_enddef(ncid);
        TEST_HANDLE_ERR(status);
        MPI_Barrier(MPI_COMM_WORLD);
        def_time = MPI_Wtime()-start_time-open_time;

/* to eliminate paging effects, do the operations once but don't time
   them */
        if (mvar_flag == 0) {
            status = ncmpi_begin_indep_data(ncid);
            TEST_HANDLE_ERR(status);
            for (i=0; i<nvars; i++){
                status = ncmpi_put_vara(ncid, varid[i],
                                        starts_list[i], count_list[i],
                                        (const void *)&(buf[i][0]),
                                        bufcount_list[i], MPI_INT);
                TEST_HANDLE_ERR(status);
            }
            status = ncmpi_end_indep_data(ncid);
            TEST_HANDLE_ERR(status);
        }
        if (mvar_flag == 1) {
            status = ncmpi_mput_vara_all(ncid, nvars, varid,
                                         starts_list, count_list,
                                         (void **)buf, bufcount_list,
                                         datatype_list);
            TEST_HANDLE_ERR(status);
        }
        if (mvar_flag == 2) {
            for (i=0; i<nvars; i++){
                status = ncmpi_iput_vara(ncid, varid[i], starts_list[i],
                                         count_list[i], (void *)&(buf[i][0]),
                                         bufcount_list[i], MPI_INT,
                                         &array_of_requests[i]);
                TEST_HANDLE_ERR(status);
                status = ncmpi_begin_indep_data(ncid);
                TEST_HANDLE_ERR(status);
                status = ncmpi_wait(ncid, 1, &array_of_requests[i],
                                    &array_of_statuses[i]);
                TEST_HANDLE_ERR(status);
                status = ncmpi_end_indep_data(ncid);
                TEST_HANDLE_ERR(status);
            }
        }
        if (mvar_flag == 3) {
            for (i=0; i<nvars; i++){
                status = ncmpi_iput_vara(ncid, varid[i], starts_list[i],
                                         count_list[i],
                                         (const void*)&(buf[i][0]),
                                         bufcount_list[i], MPI_INT,
                                         &array_of_requests[i]);
                TEST_HANDLE_ERR(status);
            }
            status = ncmpi_wait_all(ncid, nvars, array_of_requests,
                                    array_of_statuses);
            TEST_HANDLE_ERR(status);
        }
        if (mvar_flag == 4) {
            for (i=0; i<nvars; i++){
                status = ncmpi_iput_vara(ncid, varid[i], starts_list[i],
                                         count_list[i],
                                         (const void *)&(buf[i][0]),
                                         bufcount_list[i], MPI_INT,
                                         &array_of_requests[i]);
                TEST_HANDLE_ERR(status);
            }
            status = ncmpi_wait_all(ncid, nvars, array_of_requests,
                                    array_of_statuses);
            TEST_HANDLE_ERR(status);
        }
        if (mvar_flag == 5) {
            for (i=0; i<nvars; i++){
                status = ncmpi_iput_vara(ncid, varid[i], starts_list[i],
                                         count_list[i],
                                         (const void *)&(buf[i][0]),
                                         bufcount_list[i], MPI_INT,
                                         &array_of_requests[i]);
                TEST_HANDLE_ERR(status);
            }
            status = ncmpi_begin_indep_data(ncid);
            TEST_HANDLE_ERR(status);
            status = ncmpi_wait(ncid, nvars, array_of_requests,
                                array_of_statuses);
            TEST_HANDLE_ERR(status);
            status = ncmpi_end_indep_data(ncid);
            TEST_HANDLE_ERR(status);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        write_time = MPI_Wtime() - start_time - open_time - def_time;

        status = ncmpi_close(ncid);
        TEST_HANDLE_ERR(status);

        run_time = MPI_Wtime() - start_time;
        MPI_Allreduce(&open_time, &new_open_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&def_time, &new_def_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&write_time, &new_write_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&run_time, &new_run_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
    }/*end k_loop */

    if (mynod == 0) {
        int top = 0;
        for (k=0; k<(k_loop-1); k++) {
            if (new_run_tim[top]>new_run_tim[k+1])
                top = k+1;
        }

        fprintf(stderr, "mvar nvars:%d, Global array size %d x %d x %d integers, local array size: %lld x %lld x%lld\n", nvars, array_of_gsizes[0], array_of_gsizes[1], array_of_gsizes[2],lld(sizes[0]), lld(sizes[1]), lld(sizes[2]));
        fprintf(stderr, "%lldx%lldx%lld, %d: nvars:%d, loop:%d, k:%d, open_t = %f, def_t =%f, write_t = %f sec,run_t = %f sec\n", lld(sizes[0]), lld(sizes[1]), lld(sizes[2]),mvar_flag, nvars, k_loop, top, new_open_tim[top], new_def_tim[top], new_write_tim[top], new_run_tim[top]);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Info_free(&info);

    free(array_of_requests);
    free(new_open_tim);

    for (i=0; i<nvars; i++){
        free(buf[i]);
        free(starts_list[i]);
        free(count_list[i]);
    }
    free(buf);
    free(bufcount_list);
    free(varid);
    free(starts_list);
    free(count_list);
    free(pathname);

    MPI_Finalize();
    return 0;
}
