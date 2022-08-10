/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <assert.h>
#include <errno.h>

#include <mpi.h>
#include <pnetcdf.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program writes a series of 2D variables with optional data partitioning
 * patterns: block-block, *-cyclic, block-*, and *-block. The same partitioning
 * patterns are used for read after write. In both cases, nonblocking APIs are
 * used to evaluate the performance.
 *
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file.
 *
 *    % mpicc -O2 -o aggregation aggregation.c -lpnetcdf
 *
 *    % mpiexec -n 4 ./aggregation -l 5 /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *      netcdf testfile {
 *      // file format: CDF-5 (big variables)
 *      dimensions:
 *              Block_Block_Y = 10 ;
 *              Block_Block_X = 10 ;
 *              Star_Cyclic_Y = 5 ;
 *              Star_Cyclic_X = 20 ;
 *              Block_Star_Y = 20 ;
 *              Block_Star_X = 5 ;
 *              Star_Block_Y = 5 ;
 *              Star_Block_X = 20 ;
 *      variables:
 *              float block_block_var_0(Block_Block_Y, Block_Block_X) ;
 *              float star_cyclic_var_1(Star_Cyclic_Y, Star_Cyclic_X) ;
 *              float block_star_var_2(Block_Star_Y, Block_Star_X) ;
 *              float star_block_var_3(Star_Block_Y, Star_Block_X) ;
 *      data:
 *
 *       block_block_var_0 =
 *        0, 0, 0, 0, 0, 2, 2, 2, 2, 2,
 *        0, 0, 0, 0, 0, 2, 2, 2, 2, 2,
 *        0, 0, 0, 0, 0, 2, 2, 2, 2, 2,
 *        0, 0, 0, 0, 0, 2, 2, 2, 2, 2,
 *        0, 0, 0, 0, 0, 2, 2, 2, 2, 2,
 *        1, 1, 1, 1, 1, 3, 3, 3, 3, 3,
 *        1, 1, 1, 1, 1, 3, 3, 3, 3, 3,
 *        1, 1, 1, 1, 1, 3, 3, 3, 3, 3,
 *        1, 1, 1, 1, 1, 3, 3, 3, 3, 3,
 *        1, 1, 1, 1, 1, 3, 3, 3, 3, 3 ;
 *
 *       star_cyclic_var_1 =
 *        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
 *        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 ;
 *
 *       block_star_var_2 =
 *        0, 0, 0, 0, 0,
 *        0, 0, 0, 0, 0,
 *        0, 0, 0, 0, 0,
 *        0, 0, 0, 0, 0,
 *        0, 0, 0, 0, 0,
 *        1, 1, 1, 1, 1,
 *        1, 1, 1, 1, 1,
 *        1, 1, 1, 1, 1,
 *        1, 1, 1, 1, 1,
 *        1, 1, 1, 1, 1,
 *        2, 2, 2, 2, 2,
 *        2, 2, 2, 2, 2,
 *        2, 2, 2, 2, 2,
 *        2, 2, 2, 2, 2,
 *        2, 2, 2, 2, 2,
 *        3, 3, 3, 3, 3,
 *        3, 3, 3, 3, 3,
 *        3, 3, 3, 3, 3,
 *        3, 3, 3, 3, 3,
 *        3, 3, 3, 3, 3 ;
 *
 *       star_block_var_3 =
 *        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
 *        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3 ;
 *      }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ERR(e) { \
    if ((e) != NC_NOERR) { \
        printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(e)); \
        nerrs++; \
    } \
}

#define DBG_PRINT(pattern, i, j) { \
    printf("%s i=%d j=%d: start=%lld %lld count=%lld %lld\n", \
           pattern, i, j, start[0], start[1], count[0], count[1]); \
}

#define XTYPE NC_FLOAT

static int debug;

typedef struct {
    int num_records;
    int nvars;
    int block_block;
    int star_cyclic;
    int block_star;
    int star_block;
    int blocking_io;
    MPI_Offset len;
    MPI_Offset w_size;
    MPI_Offset r_size;
    MPI_Offset header_size;
    MPI_Offset header_extent;
    MPI_Info w_info_used;
    MPI_Info r_info_used;
} config;

/*----< print_info() >------------------------------------------------------*/
static
void print_info(MPI_Info *info_used)
{
    int  i, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %24s, value = %s\n",i,key,value);
    }
}

/*----< benchmark_write() >---------------------------------------------------*/
static
int benchmark_write(char       *filename,
                    config     *cfg,
                    double     *timing)  /* [6] */
{
    int i, j, k, v, n, rank, nprocs, nerrs=0, err, num_reqs, nvars;
    int ncid, cmode, *varid, *reqs, *sts, psizes[2], time_id;
    int bb_dimids[3], sc_dimids[3], bs_dimids[3], sb_dimids[3];
    double **buf;
    double start_t, end_t;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Offset gsizes[2], start[3], count[3], lenlen;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* set PnetCDF I/O hints */
    MPI_Info_create(&info);
    /* disable the header extent alignments
    MPI_Info_set(info, "nc_header_align_size", "1");    size in bytes
    */
    /* disable the fixed-size variable alignments */
    MPI_Info_set(info, "nc_var_align_size", "1");

    nvars = 0;
    if (cfg->block_block) nvars++;
    if (cfg->star_cyclic) nvars++;
    if (cfg->block_star) nvars++;
    if (cfg->star_block) nvars++;
    nvars *= cfg->nvars;

    varid = (int*) malloc(nvars * sizeof(int));

    /* initialize I/O buffer */
    lenlen = cfg->len * cfg->len;
    buf = (double**) malloc(nvars * sizeof(double*));
    for (i=0; i<nvars; i++) {
        buf[i] = (double*) malloc(lenlen * sizeof(double*));
        assert(buf[i] != NULL);
        for (j=0; j<lenlen; j++) buf[i][j] = (double)rank;
    }
    MPI_Barrier(comm);
    timing[0] = MPI_Wtime();

    /* create a new file for writing -----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR(err)
    start_t = MPI_Wtime();
    timing[1] = start_t - timing[0];
    MPI_Info_free(&info);

    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &time_id); ERR(err)
    if (cfg->block_block) {
        gsizes[0] = cfg->len * psizes[0];
        gsizes[1] = cfg->len * psizes[1];
        err = ncmpi_def_dim(ncid, "Block_Block_Y", gsizes[0], &bb_dimids[1]);
        ERR(err)
        err = ncmpi_def_dim(ncid, "Block_Block_X", gsizes[1], &bb_dimids[2]);
        ERR(err)
        bb_dimids[0] = time_id;
    }
    if (cfg->star_cyclic) {
        gsizes[0] = cfg->len;
        gsizes[1] = cfg->len * nprocs;
        err = ncmpi_def_dim(ncid, "Star_Cyclic_Y", gsizes[0], &sc_dimids[1]);
        ERR(err)
        err = ncmpi_def_dim(ncid, "Star_Cyclic_X", gsizes[1], &sc_dimids[2]);
        ERR(err)
        sc_dimids[0] = time_id;
    }
    if (cfg->block_star) {
        gsizes[0] = cfg->len * nprocs;
        gsizes[1] = cfg->len;
        err = ncmpi_def_dim(ncid, "Block_Star_Y",  gsizes[0], &bs_dimids[1]);
        ERR(err)
        err = ncmpi_def_dim(ncid, "Block_Star_X",  gsizes[1], &bs_dimids[2]);
        ERR(err)
        bs_dimids[0] = time_id;
    }
    if (cfg->star_block) {
        gsizes[0] = cfg->len;
        gsizes[1] = cfg->len * nprocs;
        err = ncmpi_def_dim(ncid, "Star_Block_Y",  gsizes[0], &sb_dimids[1]);
        ERR(err)
        err = ncmpi_def_dim(ncid, "Star_Block_X",  gsizes[1], &sb_dimids[2]);
        ERR(err)
        sb_dimids[0] = time_id;
    }

    /* define variables */
    v = num_reqs = 0;
    for (i=0; i<cfg->nvars; i++) {
        char name[32];
        if (cfg->block_block) {
            /* variables are block-block partitioned */
            sprintf(name,"block_block_var_%d",v);
            err = ncmpi_def_var(ncid, name, XTYPE, 3, bb_dimids, &varid[v++]);
            ERR(err)
            num_reqs++;
        }
        if (cfg->star_cyclic) {
            /* variables are *-cyclic partitioned */
            sprintf(name,"star_cyclic_var_%d",v);
            err = ncmpi_def_var(ncid, name, XTYPE, 3, sc_dimids, &varid[v++]);
            ERR(err)
            num_reqs += cfg->len;
        }
        if (cfg->block_star) {
            /* variables are block-* partitioned */
            sprintf(name,"block_star_var_%d",v);
            err = ncmpi_def_var(ncid, name, XTYPE, 3, bs_dimids, &varid[v++]);
            ERR(err)
            num_reqs++;
        }
        if (cfg->star_block) {
            /* variables are *-block partitioned */
            sprintf(name,"star_block_var_%d",v);
            err = ncmpi_def_var(ncid, name, XTYPE, 3, sb_dimids, &varid[v++]);
            ERR(err)
            num_reqs++;
        }
    }
    assert(v == nvars);
    reqs = (int*) malloc(num_reqs * sizeof(int));
    sts  = (int*) malloc(num_reqs * sizeof(int));

    err = ncmpi_enddef(ncid); ERR(err)
    err = ncmpi_inq_header_size(ncid, &cfg->header_size); ERR(err)
    err = ncmpi_inq_header_extent(ncid, &cfg->header_extent); ERR(err)
    end_t = MPI_Wtime();
    timing[2] = end_t - start_t;
    start_t = end_t;

    timing[3] = timing[4] = 0;
    for (n=0; n<cfg->num_records; n++) {
        start_t = MPI_Wtime();
        k = v = 0;
        start[0] = n;
        count[0] = 1;
        for (i=0; i<cfg->nvars; i++) {
            if (cfg->block_block) {
                start[1] = cfg->len * (rank % psizes[0]);
                start[2] = cfg->len * ((rank / psizes[1]) % psizes[1]);
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_put_vara_double_all(ncid, varid[v], start,
                                                    count, buf[v]);
                else
                    err = ncmpi_iput_vara_double(ncid, varid[v], start, count,
                                                 buf[v], &reqs[k++]);
                ERR(err)
                if (debug) DBG_PRINT("block-block", i, 0)
                v++;
            }
            if (cfg->star_cyclic) {
                double *ptr = buf[v];
                start[1] = 0;
                start[2] = rank;
                count[1] = cfg->len;
                count[2] = 1;
                for (j=0; j<cfg->len; j++) {
                    if (cfg->blocking_io)
                        err = ncmpi_put_vara_double_all(ncid, varid[v], start,
                                                        count, ptr);
                    else
                        err = ncmpi_iput_vara_double(ncid, varid[v], start, count,
                                                     ptr, &reqs[k++]);
                    ERR(err)
                    ptr += cfg->len;
                    start[2] += nprocs;
                    if (debug) DBG_PRINT("*-cyclic", i, j)
                }
                v++;
            }
            if (cfg->block_star) {
                start[1] = cfg->len * rank;
                start[2] = 0;
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_put_vara_double_all(ncid, varid[v], start,
                                                    count, buf[v]);
                else
                    err = ncmpi_iput_vara_double(ncid, varid[v], start, count,
                                                 buf[v], &reqs[k++]);
                ERR(err)
                if (debug) DBG_PRINT("block-*", i, 0)
                v++;
            }
            if (cfg->star_block) {
                start[1] = 0;
                start[2] = cfg->len * rank;
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_put_vara_double_all(ncid, varid[v], start,
                                                    count, buf[v]);
                else
                    err = ncmpi_iput_vara_double(ncid, varid[v], start, count,
                                                 buf[v], &reqs[k++]);
                ERR(err)
                if (debug) DBG_PRINT("*-block", i, 0)
                v++;
            }
        }
        assert(nvars == v);
        if (!cfg->blocking_io) assert(num_reqs == k);

        end_t = MPI_Wtime();
        timing[3] += end_t - start_t;

        if (!cfg->blocking_io) {
            start_t = end_t;
#ifdef USE_INDEP_MODE
            err = ncmpi_begin_indep_data(ncid);          ERR(err)
            err = ncmpi_wait(ncid, num_reqs, reqs, sts); ERR(err)
            err = ncmpi_end_indep_data(ncid);            ERR(err)
#else
            err = ncmpi_wait_all(ncid, num_reqs, reqs, sts); ERR(err)
#endif
            /* check status of all requests */
            for (i=0; i<num_reqs; i++) ERR(sts[i])

            end_t = MPI_Wtime();
            timing[4] += end_t - start_t;
        }
    }
    start_t = MPI_Wtime();

    /* get the true I/O amount committed */
    err = ncmpi_inq_put_size(ncid, &cfg->w_size); ERR(err)

    /* get all the hints used */
    err = ncmpi_inq_file_info(ncid, &cfg->w_info_used); ERR(err)

    err = ncmpi_close(ncid); ERR(err)

    end_t = MPI_Wtime();
    timing[5] = end_t - start_t;
    timing[0] = end_t - timing[0];

    free(sts);
    free(reqs);
    free(varid);
    for (i=0; i<nvars; i++) free(buf[i]);
    free(buf);

    return nerrs;
}

/*----< benchmark_read() >---------------------------------------------------*/
static
int benchmark_read(char       *filename,
                   config     *cfg,
                   double     *timing)  /* [5] */
{
    int i, j, k, v, n, rank, nprocs, nerrs=0, err, num_reqs, nvars;
    int ncid, *reqs, *sts, psizes[2];
    double **buf;
    double start_t, end_t;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Offset start[3], count[3], lenlen;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    lenlen = cfg->len * cfg->len;
    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    nvars = 0;
    if (cfg->block_block) nvars++;
    if (cfg->star_cyclic) nvars++;
    if (cfg->block_star) nvars++;
    if (cfg->star_block) nvars++;
    nvars *= cfg->nvars;

    /* allocate I/O buffer */
    buf = (double**) malloc(nvars * sizeof(double*));
    for (i=0; i<nvars; i++) {
        buf[i] = (double*) malloc(lenlen * sizeof(double*));
        assert(buf[i] != NULL);
    }
    MPI_Barrier(comm);
    timing[0] = MPI_Wtime();

    /* open file for reading -----------------------------------------*/
    err = ncmpi_open(comm, filename, NC_NOWRITE, info, &ncid); ERR(err)
    start_t = MPI_Wtime();
    timing[1] = start_t - timing[0];

    /* Note that PnetCDF read the file in chunks of size 256KB, thus the read
     * amount may be more than the file header size
     */
    err = ncmpi_inq_header_size(ncid, &cfg->header_size); ERR(err)
    err = ncmpi_inq_header_extent(ncid, &cfg->header_extent); ERR(err)

    num_reqs = 0;
    for (i=0; i<cfg->nvars; i++) {
        if (cfg->block_block)
            num_reqs++; /* complete in 1 nonblocking call */
        if (cfg->star_cyclic)
            num_reqs += cfg->len; /* complete in cfg->len nonblocking calls */
        if (cfg->block_star)
            num_reqs++; /* complete in 1 nonblocking call */
        if (cfg->star_block)
            num_reqs++; /* complete in 1 nonblocking call */
    }
    reqs = (int*) malloc(num_reqs * sizeof(int));
    sts  = (int*) malloc(num_reqs * sizeof(int));

    timing[2] = timing[3] = 0;
    for (n=0; n<cfg->num_records; n++) {
        start_t = MPI_Wtime();
        k = v = 0;
        start[0] = n;
        count[0] = 1;
        for (i=0; i<cfg->nvars; i++) {
            if (cfg->block_block) {
                start[1] = cfg->len * (rank % psizes[0]);
                start[2] = cfg->len * ((rank / psizes[1]) % psizes[1]);
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_get_vara_double_all(ncid, v, start, count, buf[v]);
                else
                    err = ncmpi_iget_vara_double(ncid, v, start, count, buf[v],
                                                 &reqs[k++]);
                ERR(err)
                v++;
            }
            if (cfg->star_cyclic) {
                double *ptr = buf[v];
                start[1] = 0;
                start[2] = rank;
                count[1] = cfg->len;
                count[2] = 1;
                for (j=0; j<cfg->len; j++) {
                    if (cfg->blocking_io)
                        err = ncmpi_get_vara_double_all(ncid, v, start, count, ptr);
                    else
                        err = ncmpi_iget_vara_double(ncid, v, start, count, ptr,
                                                     &reqs[k++]);
                    ERR(err)
                    ptr += cfg->len;
                    start[2] += nprocs;
                }
                v++;
            }
            if (cfg->block_star) {
                start[1] = cfg->len * rank;
                start[2] = 0;
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_get_vara_double_all(ncid, v, start, count, buf[v]);
                else
                    err = ncmpi_iget_vara_double(ncid, v, start, count, buf[v],
                                                 &reqs[k++]);
                ERR(err)
                v++;
            }
            if (cfg->star_block) {
                start[1] = 0;
                start[2] = cfg->len * rank;
                count[1] = cfg->len;
                count[2] = cfg->len;
                if (cfg->blocking_io)
                    err = ncmpi_get_vara_double_all(ncid, v, start, count, buf[v]);
                else
                    err = ncmpi_iget_vara_double(ncid, v, start, count, buf[v],
                                                 &reqs[k++]);
                ERR(err)
                v++;
            }
        }
        assert(nvars == v);
        if (!cfg->blocking_io) assert(num_reqs == k);

        end_t = MPI_Wtime();
        timing[2] += end_t - start_t;

        if (!cfg->blocking_io) {
            start_t = end_t;

#ifdef USE_INDEP_MODE
            err = ncmpi_begin_indep_data(ncid);          ERR(err)
            err = ncmpi_wait(ncid, num_reqs, reqs, sts); ERR(err)
            err = ncmpi_end_indep_data(ncid);            ERR(err)
#else
            err = ncmpi_wait_all(ncid, num_reqs, reqs, sts); ERR(err)
#endif
            /* check status of all requests */
            for (i=0; i<num_reqs; i++) ERR(sts[i])

            end_t = MPI_Wtime();
            timing[3] += end_t - start_t;
        }
    }
    start_t = MPI_Wtime();

    /* get the true I/O amount committed */
    err = ncmpi_inq_get_size(ncid, &cfg->r_size); ERR(err)

    /* get all the hints used */
    err = ncmpi_inq_file_info(ncid, &cfg->r_info_used); ERR(err)

    err = ncmpi_close(ncid); ERR(err)

    end_t = MPI_Wtime();
    timing[4] = end_t - start_t;
    timing[0] = end_t - timing[0];

    free(sts);
    free(reqs);
    for (i=0; i<nvars; i++) free(buf[i]);
    free(buf);

    return nerrs;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode\n"
    "       [-d] Debug mode\n"
    "       [-x] disable aggregation (using blocking APIs instead)\n"
    "       [-r]  read-only benchmark\n"
    "       [-w] write-only benchmark\n"
    "       [-b] block-block partitioning pattern\n"
    "       [-c] *-cyclic    partitioning pattern\n"
    "       [-i] block-*     partitioning pattern\n"
    "       [-j] *-block     partitioning pattern\n"
    "       [-l len]: local variable of size len x len (default 10)\n"
    "       [-n num]: number of variables each pattern (default 1)\n"
    "       [-t num]: number of time records (default 1)\n"
    "       [filename]: output netCDF file name (default ./testfile.nc)\n\n"
    " When both -r and -w are not set, write and read benchmarks are enabled\n"
    " When none of pattern options is set, all patterns are enabled\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char** argv) {
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, nprocs, verbose=1, nerrs=0, enable_read, enable_write;
    int nvars, block_block, star_cyclic, block_star, star_block, num_records;
    int blocking_io;
    double timing[11], max_t[11];
    MPI_Offset len=0, sum_w_size, sum_r_size;
    MPI_Comm comm=MPI_COMM_WORLD;
    config cfg;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    nvars        = 1;
    block_block  = 0;
    star_cyclic  = 0;
    block_star   = 0;
    star_block   = 0;
    enable_read  = 0;
    enable_write = 0;
    num_records  = 1;
    blocking_io  = 0;

    /* get command-line arguments */
    debug = 0;
    while ((i = getopt(argc, argv, "hqdbcijrwxl:n:t:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'd': debug = 1;
                      break;
            case 'b': block_block = 1;
                      break;
            case 'c': star_cyclic = 1;
                      break;
            case 'i': block_star = 1;
                      break;
            case 'j': star_block = 1;
                      break;
            case 'x': blocking_io = 1;
                      break;
            case 'r': enable_read = 1;
                      break;
            case 'w': enable_write = 1;
                      break;
            case 'l': len = atoi(optarg);
                      break;
            case 'n': nvars = atoi(optarg);
                      break;
            case 't': num_records = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    len = (len <= 0) ? 10 : len;

    if (block_block == 0 && star_cyclic == 0 && block_star == 0 &&
        star_block == 0)
        block_block = star_cyclic = block_star = star_block = 1;

    cfg.nvars       = nvars;
    cfg.block_block = block_block;
    cfg.star_cyclic = star_cyclic;
    cfg.block_star  = block_star;
    cfg.star_block  = star_block;
    cfg.len         = len;
    cfg.num_records = num_records;
    cfg.blocking_io = blocking_io;

    if (enable_read == 0 && enable_write == 0)
        enable_read = enable_write = 1;

    if (enable_write)
        nerrs += benchmark_write(filename, &cfg, timing);
    if (enable_read)
        nerrs += benchmark_read (filename, &cfg, timing+6);

    MPI_Reduce(&timing,     &max_t,     11, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&cfg.w_size, &sum_w_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
    MPI_Reduce(&cfg.r_size, &sum_r_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
    if (verbose && rank == 0) {
        double bw;
        printf("-----------------------------------------------------------\n");
        print_info(&cfg.w_info_used);
        printf("-----------------------------------------------------------\n");
        nvars = 0;
        if (cfg.block_block) {
            printf("benchmarking block-block partitioning pattern: enabled\n");
            nvars++;
        }
        if (cfg.star_cyclic) {
            printf("benchmarking *-cyclic partitioning pattern:    enabled\n");
            nvars++;
        }
        if (cfg.block_star) {
            printf("benchmarking block-* partitioning pattern:     enabled\n");
            nvars++;
        }
        if (cfg.star_block) {
            printf("benchmarking *-block partitioning pattern:     enabled\n");
            nvars++;
        }
        printf("-----------------------------------------------------------\n");
        nvars *= cfg.nvars;
        printf("Output NetCDF file name: %s\n", filename);
        printf("Output NetCDF file header size:         %lld B\n", cfg.header_size);
        printf("Output NetCDF file header extent:       %lld B\n", cfg.header_extent);
        printf("Number of MPI processes:                %d\n", nprocs);
        printf("Total number of variables:              %d\n", nvars);
        if (XTYPE == NC_FLOAT)
            printf("Data type of variables in output file:  NC_FLOAT\n");
        else if (XTYPE == NC_DOUBLE)
            printf("Data type of variables in output file:  NC_DOUBLE\n");
        printf("Data type of variables in memory:       double\n");
        printf("Local 2D variable size in each process: %lld x %lld\n",len,len);
        printf("Number of time records:                 %d\n",num_records);
        printf("-----------------------------------------------------------\n");
        if (enable_write) {
            bw = (double)sum_w_size / 1048576.0;
            printf("Total write amount        = %11lld B = %9.2f MiB = %6.2f GiB\n", sum_w_size, bw, bw/1024);
            printf("Max file open/create time = %16.4f sec\n", max_t[1]);
            printf("Max PnetCDF define   time = %16.4f sec\n", max_t[2]);
            if (cfg.blocking_io)
                printf("Max   blocking write time = %16.4f sec\n", max_t[3]);
            else {
                printf("Max nonblocking post time = %16.4f sec\n", max_t[3]);
                printf("Max nonblocking wait time = %16.4f sec\n", max_t[4]);
            }
            printf("Max file close       time = %16.4f sec\n", max_t[5]);
            printf("Max open-to-close    time = %16.4f sec\n", max_t[0]);
            printf("Write bandwidth           = %14.2f MiB/s = %9.2f GiB/s\n", bw/max_t[0], bw/1024.0/max_t[0]);
            printf("-------------------------------------------------------\n");
        }
        if (enable_read) {
            bw = (double)sum_r_size / 1048576.0;
            printf("Total read  amount        = %11lld B = %9.2f MiB = %6.2f GiB\n", sum_r_size, bw, bw/1024);
            printf("Max file open/create time = %16.4f sec\n", max_t[7]);
            if (cfg.blocking_io)
                printf("Max   blocking  read time = %16.4f sec\n", max_t[7]);
            else {
                printf("Max nonblocking post time = %16.4f sec\n", max_t[8]);
                printf("Max nonblocking wait time = %16.4f sec\n", max_t[9]);
            }
            printf("Max file close       time = %16.4f sec\n", max_t[10]);
            printf("Max open-to-close    time = %16.4f sec\n", max_t[6]);
            printf("Read  bandwidth           = %14.2f MiB/s = %9.2f GiB/s\n", bw/max_t[6], bw/1024.0/max_t[6]);
            printf("-------------------------------------------------------\n");
        }
    }
    if (enable_write) MPI_Info_free(&cfg.w_info_used);
    if (enable_read)  MPI_Info_free(&cfg.r_info_used);

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
    /* report the PnetCDF internal heap memory allocation high water mark */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
            printf("Max heap memory allocated by PnetCDF internally is %.2f MiB\n",
                   (float)sum_size/1048576);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

