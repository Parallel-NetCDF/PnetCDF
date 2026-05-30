/*********************************************************************
 *
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  This examples reads the variable FLDS from an output file generated from a
 *  production run of E3SM Land Model, then compress it using ZLIB and writes
 *  to a new file. This program is designed to demonstrate the usage of the
 *  chunking-compression feature of PnetCDF. Note the input file can also be a
 *  chunked-compressed PnetCDF file.
 *
 * The FLDS input file has the following metadata.
 *
 * // file format: CDF-5 (big variables)
 *  dimensions:
 *      x = 7814 ;
 *      y = 8075 ;
 *      time = UNLIMITED ; // (248 currently)
 *  variables:
 *      float x(x) ;
 *      float y(y) ;
 *      float time(time) ;
 *      float FLDS(time, y, x) ;
 *          FLDS:long_name = "incident longwave radiation" ;
 *          FLDS:units = "W/m**2" ;
 *  }
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

static int verbose;

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, \
               ncmpi_strerror(err)); \
        nerrs++; \
        goto err_out; \
    } \
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h | -q | -t | -c] [-k format] -i in_file -o out_file]n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 5 for CDF-5\n"
    "       [-t]: data partitioning along time dimension (default: no)\n"
    "       [-c]: use cyclic partitioning pattern, only relevant when -t is used (default: block)\n"
    "       -i filename:  input netCDF file name\n"
    "       -o filename: output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
pnetcdf_io(MPI_Comm    comm,
           const char *in_path,
           const char *out_path,
           int         cmode,
           int         div_time,
           int         parti)
{
    int i, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[3], ntimes, chunk_dim[3];
    float *buf = NULL, *buf_ptr;
    double timing[2], max_t[2];
    MPI_Offset tlen, ylen, xlen, start[3], count[3], amnt[2], sum_amnt[2];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");
    MPI_Info_set(info, "nc_chunk_default_filter", "zlib");

    /* chunking is supported only when MPI-IO driver is used */
    MPI_Info_set(info, "nc_pncio", "disable");

    /* open the input file file */
    err = ncmpi_open(comm, in_path, NC_NOWRITE, info, &ncid); ERR
    err = ncmpi_inq_dimid(ncid, "time", &dimid[0]); ERR
    err = ncmpi_inq_dimid(ncid, "y",    &dimid[1]); ERR
    err = ncmpi_inq_dimid(ncid, "x",    &dimid[2]); ERR

    err = ncmpi_inq_dimlen(ncid, dimid[0], &tlen); ERR
    err = ncmpi_inq_dimlen(ncid, dimid[1], &ylen); ERR
    err = ncmpi_inq_dimlen(ncid, dimid[2], &xlen); ERR

    err = ncmpi_inq_varid(ncid, "FLDS", &varid); ERR

    if (div_time) {
        /* partition along time dimension */
        chunk_dim[0] = 1;
        chunk_dim[1] = ylen;
        chunk_dim[2] = xlen;

        ntimes = tlen / nprocs;
        if (rank < tlen % nprocs)
            ntimes++;

        if (parti) { /* block partitioning */
            start[0] = (tlen / nprocs) * rank;
            if (rank < tlen % nprocs)
                start[0] += rank;
            else
                start[0] += tlen % nprocs;
        }
        else
            start[0] = rank; /* cyclic partitioning */

        start[1] = 0;
        start[2] = 0;

        count[0] = 1;
        count[1] = ylen;
        count[2] = xlen;
    }
    else {
        /* checkerboard partitioning on every time step */
        int psize[2], yrank, xrank;

        chunk_dim[0] = 1;
        chunk_dim[1] = 1010;
        chunk_dim[2] = 977;

        /* Creates a division of processors in a Cartesian grid */
        psize[0] = psize[1] = 0;
        MPI_Dims_create(nprocs, 2, psize);

        yrank = rank / psize[1];
        xrank = rank % psize[1];

        if (verbose) {
            if (rank == 0) printf("psize %d %d\n", psize[0],psize[1]);
            printf("%2d: yrank %d xrank %d\n", rank, yrank, xrank);
        }

        ntimes = tlen;

        start[0] = 0;
        count[0] = 1;

        count[1] = ylen / psize[0];
        start[1] = count[1] * yrank;
        if (yrank < ylen % psize[0]) {
            start[1] += yrank;
            count[1]++;
        }
        else
            start[1] += ylen % psize[0];

        count[2] = xlen / psize[1];
        start[2] = count[2] * xrank;
        if (xrank < xlen % psize[1]) {
            start[2] += xrank;
            count[2]++;
        }
        else
            start[2] += xlen % psize[1];
    }

    if (verbose) {
        printf("%2d: ntimes %d start %lld %lld %lld count %lld %lld %lld end %lld %lld\n",
               rank,ntimes,start[0],start[1],start[2],count[0],count[1],count[2],
               start[1]+count[1],start[2]+count[2]);
        fflush(stdout);
    }

    /* allocate read buffer */
    buf = (float*) malloc(sizeof(float) * ntimes*count[1]*count[2]);

    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();

    buf_ptr = buf;
    for (i=0; i<ntimes; i++) {

        /* use PnetCDF nonblocking get API */
        err = ncmpi_iget_vara_float(ncid, varid, start, count, buf_ptr, NULL);
        ERR

        buf_ptr += count[1]*count[2];

        if (div_time) {
            if (parti)
                start[0]++; /* block partitioning */
            else
                start[0] += nprocs; /* cyclic partitioning */
        }
        else
            start[0]++;
    }

    /* collectively wait for reads to complete */
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR
    timing[0] = MPI_Wtime() - timing[0];

    err = ncmpi_inq_get_size(ncid, &amnt[0]); ERR

    err = ncmpi_close(ncid); ERR

    if (verbose && rank == 0)
        printf("---- read completed -----------------------------\n");

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, out_path, cmode, info, &ncid); ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "y",    ylen,         &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "x",    xlen,         &dimid[2]); ERR

    /* define variable */
    err = ncmpi_def_var(ncid, "FLDS", NC_FLOAT, 3, dimid, &varid); ERR

    /* set chunking (1st dimension should always be 1 for record variable) */
    err = ncmpi_var_set_chunk(ncid, varid, chunk_dim); ERR

    /* use default filter */
    err = ncmpi_var_set_filter(ncid, varid, NC_FILTER_DEFLATE); ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    if (div_time) {
        if (parti) { /* block partitioning */
            start[0] = (tlen / nprocs) * rank;
            if (rank < tlen % nprocs)
                start[0] += rank;
            else
                start[0] += tlen % nprocs;
        }
        else
            start[0] = rank; /* cyclic partitioning */
    }
    else
        start[0] = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    timing[1] = MPI_Wtime();

    buf_ptr = buf;
    for (i=0; i<ntimes; i++) {
        /* use PnetCDF nonblocking put API */
        err = ncmpi_iput_vara_float(ncid, varid, start, count, buf_ptr, NULL);
        ERR

        buf_ptr += count[1]*count[2];

        if (div_time) {
            if (parti)
                start[0]++; /* block partitioning */
            else
                start[0] += nprocs; /* cyclic partitioning */
        }
        else
            start[0]++;
    }

    /* collectively wait for writes to complete */
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR
    timing[1] = MPI_Wtime() - timing[1];

    err = ncmpi_inq_put_size(ncid, &amnt[1]); ERR

    /* close the file */
    err = ncmpi_close(ncid); ERR

    MPI_Reduce(amnt, sum_amnt, 2, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(timing, max_t, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        MPI_Offset vsize = sizeof(float) * tlen * ylen * xlen;
        printf("Data partitioning along: %s\n", (div_time)?"time dimension":"Y-X dimensions");
        if (div_time)
            printf("Data partitioning pattern: %s\n", (parti)?"CYCLIC":"BLOCK");
        printf("Dimension time of size:  %10lld\n", tlen);
        printf("Dimension y of size:     %10lld\n", ylen);
        printf("Dimension x of size:     %10lld\n", xlen);
        printf("Global variable size:    %lld bytes\n", vsize);
        printf("Total read  data amount: %lld bytes\n", sum_amnt[0]);
        printf("Total write data amount: %lld bytes\n", sum_amnt[1]);
        printf("Max  read timing:        %10.2f seconds\n", max_t[0]);
        printf("Read   bandwidth:        %10.2f MiB/s\n", (double)sum_amnt[0]/1048576.0/max_t[0]);
        printf("Max write timing:        %10.2f seconds\n", max_t[1]);
        printf("Write  bandwidth:        %10.2f MiB/s\n\n", (double)sum_amnt[1]/1048576.0/max_t[1]);
    }

err_out:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);
    if (buf != NULL) free(buf);

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char *in_path=NULL, *out_path=NULL;
    int i, rank, kind=0, cmode=0, nerrs=0, parti, div_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;
    div_time = 0;
    parti = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqtck:i:o:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 't': div_time = 1; /* divide along time dimension */
                      break;
            case 'c': parti = 1; /* cyclic partitioning */
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'i': in_path = strdup(optarg);
                      break;
            case 'o': out_path = strdup(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (in_path == NULL || out_path == NULL) {
        if (rank == 0) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }

    if (verbose && rank == 0) printf("%s: test compression using ELM's FLDS data\n",__FILE__);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET; break;
        case(5): cmode = NC_64BIT_DATA;   break;
        default: cmode = 0;
    }

    nerrs += pnetcdf_io(MPI_COMM_WORLD, in_path, out_path, cmode, div_time, parti);

    free(in_path);
    free(out_path);

    MPI_Finalize();
    return (nerrs > 0);
}

