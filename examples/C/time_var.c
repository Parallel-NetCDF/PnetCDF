/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to create and write record and fixed-size variables.
 * It first defines a record 2D variable of size time * global_nx where
 *    time is a expandable dimension and
 *    global_nx == (NX * number of MPI processes).
 * The data partitioning pattern is a column-wise partitioning across all
 * processes. Each process writes a subarray of size 1 * nx.
 * It then defines a fixed-size 2D variable of size global_ny * global_nx where
 *    global_ny == NY and
 *    global_nx == (NX * number of MPI processes).
 * The data partitioning pattern is a column-wise partitioning across all
 * processes. Each process writes a subarray of size ny * nx.
 *
 *    To compile:
 *        mpicc -O2 time_var.c -o time_var -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output NetCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./time_var /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *            time = UNLIMITED ; // (2 currently)
 *            Y = 4 ;
 *            X = 12 ;
 *    variables:
 *            float rec_var(time, X) ;
 *            float fix_var(Y, X) ;
 *
 *    // global attributes:
 *                    :history = "Mon Aug 13 21:27:48 2018" ;
 *    data:
 *
 *     rec_var =
 *      100, 100, 100, 101, 101, 101, 102, 102, 102, 103, 103, 103,
 *      200, 200, 200, 201, 201, 201, 202, 202, 202, 203, 203, 203 ;
 *
 *     fix_var =
 *      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
 *      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
 *      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
 *      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 4
#define NX 3

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 3 for NetCDF4,\n"
    "                                4 for NetCDF4 classic model, 5 for CDF-5\n"
    "       [filename] output netCDF file name\n";
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

/*----< pnetcdf_write() >----------------------------------------------------*/
static int
pnetcdf_write(MPI_Comm comm, char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, dimids[2], dim_t, dim_y, dim_x, rec_var, fix_var;
    double buf[NY][NX];
    char str_att[128];
    MPI_Offset global_ny, global_nx;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));
        return 1;
    }

    /* set the global dimensions NY and (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * nprocs;

    /* add a global attribute: a time stamp at rank 0 */
    time_t ltime = time(NULL); /* get the current calendar time */
    asctime_r(localtime(&ltime), str_att);
    sprintf(str_att, "Mon Aug 13 21:27:48 2018");

    /* make sure the time string are consistent among all processes */
    MPI_Bcast(str_att, sizeof(str_att), MPI_CHAR, 0, comm);

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             &str_att[0]); ERR

    /* define dimensions time, Y, and X */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dim_t); ERR
    err = ncmpi_def_dim(ncid, "Y",    global_ny,    &dim_y); ERR
    err = ncmpi_def_dim(ncid, "X",    global_nx,    &dim_x); ERR

    /* define a 2D record variable of float type */
    dimids[0] = dim_t;
    dimids[1] = dim_x;
    err = ncmpi_def_var(ncid, "rec_var", NC_FLOAT, 2, dimids, &rec_var); ERR

    /* define a 2D fixed-size variable of float type */
    dimids[0] = dim_y;
    dimids[1] = dim_x;
    err = ncmpi_def_var(ncid, "fix_var", NC_FLOAT, 2, dimids, &fix_var); ERR

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* write to record variable (the variable with expandable dimension) */
    for (j=0; j<NX; j++) buf[0][j] = 1.0 * rank + 100.0;

    start[0] = 0;
    start[1] = NX * rank;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_double_all(ncid, rec_var, start, count, &buf[0][0]);ERR

    /* write to fixed-size variable (the variable with no expandable dimension) */
    for (i=0; i<NY; i++) for (j=0; j<NX; j++)
        buf[i][j] = 1.0 * rank;

    start[0] = 0;
    start[1] = NX * rank;
    count[0] = NY;
    count[1] = NX;
    err = ncmpi_put_vara_double_all(ncid, fix_var, start, count, &buf[0][0]);ERR

    /* write a new record: 2nd record */
    for (j=0; j<NX; j++) buf[1][j] = 1.0 * rank + 200.0;

    start[0] = 1;  /* 2nd record */
    start[1] = NX * rank;
    count[0] = 1;
    count[1] = NX;
    err = ncmpi_put_vara_double_all(ncid, rec_var, start, count, &buf[1][0]);ERR

    err = ncmpi_close(ncid); ERR

    return nerrs;
}

/*----< pnetcdf_read() >-----------------------------------------------------*/
static int
pnetcdf_read(MPI_Comm comm, char *filename)
{
    int i, rank, nprocs, err, nerrs=0, local_nx;
    int ncid, dim_t, dim_y, dim_x, rec_var, fix_var;
    double *buf;
    char *str_att;
    MPI_Offset time_len, global_ny, global_nx, str_len;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* open file for reading ----------------------------------------*/
    err = ncmpi_open(comm, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));
        return 1;
    }

    /* ncmpi_open automatically enters data mode */

    /* inquire dimensions */
    err = ncmpi_inq_dimid(ncid, "time", &dim_t); ERR
    err = ncmpi_inq_dimid(ncid, "Y",    &dim_y); ERR
    err = ncmpi_inq_dimid(ncid, "X",    &dim_x); ERR
    err = ncmpi_inq_dimlen(ncid, dim_t, &time_len); ERR
    err = ncmpi_inq_dimlen(ncid, dim_y, &global_ny); ERR
    err = ncmpi_inq_dimlen(ncid, dim_x, &global_nx); ERR

    local_nx = (int)global_nx / nprocs;

    /* get global attribute history */
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, "history", &str_len); ERR
    str_att = (char*) malloc(str_len);
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, "history", str_att); ERR

    /* inquire variable IDs */
    err = ncmpi_inq_varid(ncid, "rec_var", &rec_var); ERR
    err = ncmpi_inq_varid(ncid, "fix_var", &fix_var); ERR

    /* allocate read buffer */
    buf = (double*) malloc(global_ny * global_nx * sizeof(double));
    for (i=0; i<global_ny*global_nx; i++) buf[i] = -1.0; /* initialize buffer */

    /* read 1st record of the record variable, rec_var */
    start[0] = 0;
    start[1] = (MPI_Offset)local_nx * rank;
    count[0] = 1;
    count[1] = local_nx;
    err = ncmpi_get_vara_double_all(ncid, rec_var, start, count, buf); ERR

    /* check read contents */
    for (i=0; i<local_nx; i++) {
        double expect=1.0 * rank + 100.0;
        if (buf[i] != expect) {
            printf("Read error at line %d: buf[%d] expect %f but got %f\n",
                   __LINE__, i, expect, buf[i]);
            nerrs++;
            break;
        }
    }
    for (i=0; i<global_nx; i++) buf[i] = -1.0; /* reset buffer */

    /* read fixed-size variable, fix_var */
    start[0] = 0;
    start[1] = (MPI_Offset)local_nx * rank;
    count[0] = global_ny;
    count[1] = local_nx;
    err = ncmpi_get_vara_double_all(ncid, fix_var, start, count, buf); ERR

    /* check read contents */
    for (i=0; i<global_ny*local_nx; i++) {
        double expect=1.0*rank;
        if (buf[i] != expect) {
            printf("Read error at line %d: buf[%d][%d] expect %f but got %f\n",
                   __LINE__, i/local_nx, i%local_nx, expect, buf[i]);
            nerrs++;
            break;
        }
    }
    for (i=0; i<global_ny*global_nx; i++) buf[i] = -1.0; /* reset buffer */

    /* read 2nd record of variable rec_var */
    start[0] = 1;  /* 2nd record */
    start[1] = (MPI_Offset)local_nx * rank;
    count[0] = 1;
    count[1] = local_nx;
    err = ncmpi_get_vara_double_all(ncid, rec_var, start, count, buf); ERR

    /* check read contents */
    for (i=0; i<local_nx; i++) {
        double expect=1.0 * rank + 200.0;
        if (buf[i] != expect) {
            printf("Read error at line %d: buf[%d] expect %f but got %f\n",
                   __LINE__, i, expect, buf[i]);
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); ERR

    free(str_att);
    free(buf);

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, kind=0, cmode=0, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("%s: example of using put_vara APIs\n",__FILE__);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET;             break;
        case(3): cmode = NC_NETCDF4;                  break;
        case(4): cmode = NC_NETCDF4|NC_CLASSIC_MODEL; break;
        case(5): cmode = NC_64BIT_DATA;               break;
        default: cmode = 0; /* default format is CDF-1 */
    }

#ifndef PNETCDF_DRIVER_NETCDF4
    /* netcdf4 driver is not enabled, skip */
    if (kind == 3 || kind == 4) {
        MPI_Finalize();
        return 0;
    }
#endif
    nerrs += pnetcdf_write(MPI_COMM_WORLD, filename, cmode);
    if (nerrs == 0)
        nerrs += pnetcdf_read(MPI_COMM_WORLD, filename);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

