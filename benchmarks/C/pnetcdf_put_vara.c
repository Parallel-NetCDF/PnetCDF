/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program writes a series of 3D variables with 2D block-block partitioning
 * pattern. Each variable is a record variable. The number of variables,
 * variable size, number of time records, and the NetCDF file format can be
 * customized through command-line options. In addition, option '-i', if set,
 * PnetCDF nonblocking APIs will be used to write to the file.
 *
 * To compile:
 *    % mpicc -O2 pnetcdf_put_vara.c -o pnetcdf_put_vara -lpnetcdf
 *
 * Example run command and outputs from running ncmpidump on the output netCDF
 * file produced by this example program:
 *
 *    % mpiexec -n 4 ./pnetcdf_put_vara -n 4 -t 10 -l 1024 -k 5 -i
 *    -----------------------------------------------------------
 *    Output NetCDF file name:   testfile.nc
 *    Output NetCDF file kind:   CDF-5
 *    Number of MPI processes:   4
 *    Total number of variables: 4
 *    Data type of variables:    NC_FLOAT
 *    Each 2D variable size:     2048 x 2048
 *    Number of time records:    10
 *    If using nonblocking APIs: Yes
 *    Total write amount:        671088640 B
 *                               640.00 MiB
 *                               0.62 GiB
 *    Max variable write time:   0.5402 sec
 *    Write bandwidth:           1184.69 MiB/s
 *                               1.16 GiB/s
 *    Max open-to-close time:    0.5521 sec
 *    Write bandwidth:           1159.13 MiB/s
 *                               1.13 GiB/s
 *    -----------------------------------------------------------
 *
 *    % ncmpidump -h testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *    	time = UNLIMITED ; // (10 currently)
 *    	Y = 2048 ;
 *    	X = 2048 ;
 *    variables:
 *    	float var_0(time, Y, X) ;
 *    		var_0:str_att = "some text attribute 0 type text." ;
 *    		var_0:float_att = 0.f, 1.f, 2.f, 3.f ;
 *    		var_0:short_att = 1234s ;
 *    	float var_1(time, Y, X) ;
 *    		var_1:str_att = "some text attribute 1 type text." ;
 *    		var_1:float_att = 0.f, 1.f, 2.f, 3.f ;
 *    		var_1:short_att = 1234s ;
 *    	float var_2(time, Y, X) ;
 *    		var_2:str_att = "some text attribute 2 type text." ;
 *    		var_2:float_att = 0.f, 1.f, 2.f, 3.f ;
 *    		var_2:short_att = 1234s ;
 *    	float var_3(time, Y, X) ;
 *    		var_3:str_att = "some text attribute 3 type text." ;
 *    		var_3:float_att = 0.f, 1.f, 2.f, 3.f ;
 *    		var_3:short_att = 1234s ;
 *
 *    // global attributes:
 *    		:history = "Thu Aug 29 14:29:35 2024" ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

static int verbose;

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err)); \
        nerrs++; \
    } \
}

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
pnetcdf_io(char *filename,
           int nonblocking,
           int kind,
           int len,
           int nvars,
           int ntimes)
{
    char name[128], str_att[128], *file_kind;
    int i, j, rank, nprocs, err, nerrs=0;
    int psizes[2], cmode, ncid, *varid, dimid[3];
    float **buf;
    double bw, timing[2], max_t[2];
    size_t w_size;
    MPI_Offset start[3], count[3];
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    /* each local array is of size (psizes[0] * len) x (psizes[1] * len) */
    buf = (float**) malloc(sizeof(float*) * nvars);
    for (i=0; i<nvars; i++) {
        buf[i] = (float*) malloc(sizeof(float) * len * len);
        for (j=0; j<len*len; j++)
            buf[i][j] = rank + i * len + j;
    }
    varid = (int*) malloc(sizeof(int) * nvars);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET;
                 file_kind="Output NetCDF file kind:   CDF-2";
                 break;
        case(5): cmode = NC_64BIT_DATA;
                 file_kind="Output NetCDF file kind:   CDF-5";
                 break;
        default: cmode = 0;
                 file_kind="Output NetCDF file kind:   CDF-0";
    }

    /* start the timer */
    if (nprocs > 1) MPI_Barrier(comm);
    timing[0] = MPI_Wtime();

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); ERR

    /* add a global attribute */
    sprintf(str_att, "Thu Aug 29 14:29:35 2024");
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             str_att); ERR

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y", psizes[0] * len, &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X", psizes[1] * len, &dimid[2]); ERR

    /* define 2D variables of float type */
    for (i=0; i<nvars; i++) {
        sprintf(name, "var_%d", i);
        err = ncmpi_def_var(ncid, name, NC_FLOAT, 3, dimid, &varid[i]);
        ERR
    }

    /* add attributes to each variable */
    for (i=0; i<nvars; i++) {
        short short_att=1234;
        float float_att[4];

        sprintf(str_att, "some text attribute %d type text.", i);
        err = ncmpi_put_att_text(ncid, varid[i], "str_att", strlen(str_att),
                                 str_att); ERR

        for (j=0; j<4; j++) float_att[j] = j;
        err = ncmpi_put_att_float(ncid, varid[i], "float_att", NC_FLOAT, 4,
                                  float_att); ERR

        err = ncmpi_put_att_short(ncid, varid[i], "short_att", NC_SHORT, 1,
                                  &short_att); ERR
    }

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* set up the access pattern */
    start[0] = 0;
    start[1] = len * (rank / psizes[1]);
    start[2] = len * (rank % psizes[1]);
    count[0] = 1;
    count[1] = len;
    count[2] = len;

    timing[1] = MPI_Wtime();
    for (j=0; j<ntimes; j++) {
        start[0] = j;
        for (i=0; i<nvars; i++) {
            if (nonblocking)
                err = ncmpi_iput_vara_float(ncid, varid[i], start,
                                            count, buf[i], NULL);
            else
                err = ncmpi_put_vara_float_all(ncid, varid[i], start,
                                               count, buf[i]);
            ERR
        }
    }
    if (nonblocking) {
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
        ERR
    }
    timing[1] = MPI_Wtime() - timing[1];

    err = ncmpi_close(ncid); ERR
    timing[0] = MPI_Wtime() - timing[0];

    free(varid);
    for (i=0; i<nvars; i++) free(buf[i]);
    free(buf);

    if (!verbose) return nerrs;

    MPI_Reduce(timing, max_t, 2, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank > 0) return nerrs;

    printf("-----------------------------------------------------------\n");
    printf("Output NetCDF file name:   %s\n", filename);
    printf("%s\n", file_kind);
    printf("Number of MPI processes:   %d\n", nprocs);
    printf("Total number of variables: %d\n", nvars);
    printf("Data type of variables:    NC_FLOAT\n");
    printf("Each 2D variable size:     %d x %d\n",psizes[0]*len,psizes[1]*len);
    printf("Number of time records:    %d\n",ntimes);
    if (nonblocking)
        printf("If using nonblocking APIs: Yes\n");
    else
        printf("If using nonblocking APIs: NO\n");
    w_size = sizeof(float) * nprocs * len * len * nvars * ntimes;
    printf("Total write amount:        %zd B\n", w_size);
    printf("                           %.2f MiB\n", (float)w_size/1048576);
    printf("                           %.2f GiB\n", (float)w_size/1073741824);
    bw = (double)w_size / 1048576;
    printf("Max variable write time:   %.4f sec\n", max_t[1]);
    printf("Write bandwidth:           %.2f MiB/s\n", bw/max_t[1]);
    printf("                           %.2f GiB/s\n", bw/1024.0/max_t[1]);
    printf("Max open-to-close time:    %.4f sec\n", max_t[0]);
    printf("Write bandwidth:           %.2f MiB/s\n", bw/max_t[0]);
    printf("                           %.2f GiB/s\n", bw/1024.0/max_t[0]);
    printf("-----------------------------------------------------------\n");

    return nerrs;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h | -q | -i | -k num |-l num | -n num | -t num] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-i] Use PnetCDF nonblocking APIs\n"
    "       [-k format] file format: 1 for CDF-1,\n"
    "                                2 for CDF-2,\n"
    "                                5 for CDF-5 (default)\n"
    "       [-l] X and Y dimension lengths (default 32)\n"
    "       [-n] number of variables (default 4)\n"
    "       [-t] number of time iterations (default 1)\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[512];
    int i, rank, kind, nerrs=0;
    int nonblocking, nvars, len, ntimes;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;
    nonblocking = 0;
    kind = 5;
    nvars = 4;
    len = 32;
    ntimes = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqik:l:n:t:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'i': nonblocking = 1;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'l': len = atoi(optarg);
                      break;
            case 'n': nvars = atoi(optarg);
                      break;
            case 't': ntimes = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 512, "%s", argv[optind]);

    MPI_Bcast(filename, 512, MPI_CHAR, 0, MPI_COMM_WORLD);

    nerrs += pnetcdf_io(filename, nonblocking, kind, len, nvars, ntimes);

    MPI_Finalize();
    return (nerrs > 0);
}

