/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This sequential NetCDF-C program writes a series of 3D variables. Each
 * variable is a record variable. The number of variables, variable size,
 * number of time records, and the NetCDF file format can be customized through
 * command-line options.
 *
 * To compile:
 *    % gcc netcdf_put_vara.c -o netcdf_put_vara -I $NETCDF_DIR/include -L $NETCDF_DIR/lib -lnetcdf
 *
 * Example run command and outputs from running ncdump on the output netCDF
 * file produced by this example program:

 *    % ./netcdf_put_vara -n 4 -t 10 -l 1024 -k 5
 *    -----------------------------------------------------------
 *    Output NetCDF file name:   testfile.nc
 *    Output NetCDF file kind:   CDF-5
 *    Total number of variables: 4
 *    Data type of variables:    NC_FLOAT
 *    Each 2D variable size:     1024 x 1024
 *    Number of time records:    10
 *    Total write amount:        167772160 B
 *                               160.00 MiB
 *                               0.16 GiB
 *    variable write time:       0.3411 sec
 *    Write bandwidth:           469.10 MiB/s
 *                               0.46 GiB/s
 *    open-to-close time:        0.9867 sec
 *    Write bandwidth:           162.15 MiB/s
 *                               0.16 GiB/s
 *    -----------------------------------------------------------
 *
 *    % ncdump -h testfile.nc
 *    netcdf testfile {
 *    dimensions:
 *            time = UNLIMITED ; // (10 currently)
 *            Y = 1024 ;
 *            X = 1024 ;
 *    variables:
 *            float var_0(time, Y, X) ;
 *                    var_0:str_att = "some text attribute 0 type text." ;
 *                    var_0:float_att = 0.f, 1.f, 2.f, 3.f ;
 *                    var_0:short_att = 1234s ;
 *            float var_1(time, Y, X) ;
 *                    var_1:str_att = "some text attribute 1 type text." ;
 *                    var_1:float_att = 0.f, 1.f, 2.f, 3.f ;
 *                    var_1:short_att = 1234s ;
 *            float var_2(time, Y, X) ;
 *                    var_2:str_att = "some text attribute 2 type text." ;
 *                    var_2:float_att = 0.f, 1.f, 2.f, 3.f ;
 *                    var_2:short_att = 1234s ;
 *            float var_3(time, Y, X) ;
 *                    var_3:str_att = "some text attribute 3 type text." ;
 *                    var_3:float_att = 0.f, 1.f, 2.f, 3.f ;
 *                    var_3:short_att = 1234s ;
 *
 *    // global attributes:
 *                    :history = "Thu Aug 29 14:29:35 2024" ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <sys/time.h>
#include <netcdf.h>

static int verbose;

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, nc_strerror(err)); \
        nerrs++; \
    } \
}

static double wtime(void)
{
    double          now_time;
    struct timeval  etstart;
    struct timezone tzp;

    if (gettimeofday(&etstart, &tzp) == -1)
        perror("Error: calling gettimeofday() not successful.\n");

    now_time = ((double)etstart.tv_sec) +              /* in seconds */
               ((double)etstart.tv_usec) / 1000000.0;  /* in microseconds */
    return now_time;
}

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
netcdf_io(char *filename,
          int kind,
          int len,
          int nvars,
          int ntimes)
{
    char name[128], str_att[128], *file_kind;
    int i, j, err, nerrs=0;
    int cmode, ncid, *varid, dimid[3];
    float **buf;
    double bw, timing[2];
    size_t w_size;
    size_t start[3], count[3];

    /* each local array is of size len x len */
    buf = (float**) malloc(sizeof(float*) * nvars);
    for (i=0; i<nvars; i++) {
        buf[i] = (float*) malloc(sizeof(float) * len * len);
        for (j=0; j<len*len; j++)
            buf[i][j] = i * len + j;
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
    timing[0] = wtime();

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = nc_create(filename, cmode, &ncid); ERR

    /* add a global attribute */
    sprintf(str_att, "Thu Aug 29 14:29:35 2024");
    err = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             &str_att[0]); ERR

    /* define dimensions x and y */
    err = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = nc_def_dim(ncid, "Y", len, &dimid[1]); ERR
    err = nc_def_dim(ncid, "X", len, &dimid[2]); ERR

    /* define 2D variables of float type */
    for (i=0; i<nvars; i++) {
        sprintf(name, "var_%d", i);
        err = nc_def_var(ncid, name, NC_FLOAT, 3, dimid, &varid[i]);
        ERR
    }

    /* add attributes to each variable */
    for (i=0; i<nvars; i++) {
        short short_att=1234;
        float float_att[4];

        sprintf(str_att, "some text attribute %d type text.", i);
        err = nc_put_att_text(ncid, varid[i], "str_att", strlen(str_att),
                                 str_att); ERR

        for (j=0; j<4; j++) float_att[j] = j;
        err = nc_put_att_float(ncid, varid[i], "float_att", NC_FLOAT, 4,
                                  float_att); ERR

        err = nc_put_att_short(ncid, varid[i], "short_att", NC_SHORT, 1,
                                  &short_att); ERR
    }

    /* exit define mode */
    err = nc_enddef(ncid); ERR

    /* set up the access pattern */
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = len;
    count[2] = len;

    timing[1] = wtime();
    for (j=0; j<ntimes; j++) {
        start[0] = j;
        for (i=0; i<nvars; i++) {
            err = nc_put_vara_float(ncid, varid[i], start, count, buf[i]);
            ERR
        }
    }
    timing[1] = wtime() - timing[1];

    err = nc_close(ncid); ERR
    timing[0] = wtime() - timing[0];

    free(varid);
    for (i=0; i<nvars; i++) free(buf[i]);
    free(buf);

    if (!verbose) return nerrs;

    printf("-----------------------------------------------------------\n");
    printf("Output NetCDF file name:   %s\n", filename);
    printf("%s\n", file_kind);
    printf("Total number of variables: %d\n", nvars);
    printf("Data type of variables:    NC_FLOAT\n");
    printf("Each 2D variable size:     %d x %d\n",len,len);
    printf("Number of time records:    %d\n",ntimes);
    w_size = sizeof(float) * len * len * nvars * ntimes;
    printf("Total write amount:        %zd B\n", w_size);
    printf("                           %.2f MiB\n", (float)w_size/1048576);
    printf("                           %.2f GiB\n", (float)w_size/1073741824);
    bw = (double)w_size / 1048576;
    printf("variable write time:       %.4f sec\n", timing[1]);
    printf("Write bandwidth:           %.2f MiB/s\n", bw/timing[1]);
    printf("                           %.2f GiB/s\n", bw/1024.0/timing[1]);
    printf("open-to-close time:        %.4f sec\n", timing[0]);
    printf("Write bandwidth:           %.2f MiB/s\n", bw/timing[0]);
    printf("                           %.2f GiB/s\n", bw/1024.0/timing[0]);
    printf("-----------------------------------------------------------\n");

    return nerrs;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h | -q | -k num |-l num | -n num | -t num] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
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
    int i, kind, nerrs=0;
    int nvars, len, ntimes;

    verbose = 1;
    kind = 5;
    nvars = 4;
    len = 32;
    ntimes = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:l:n:t:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
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
            default:  usage(argv[0]);
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 512, "%s", argv[optind]);

    nerrs += netcdf_io(filename, kind, len, nvars, ntimes);

    return (nerrs > 0);
}

