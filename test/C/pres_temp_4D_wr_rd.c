/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
   This is an example program which writes some 4D pressure and
   temperatures. It is intended to illustrate the use of the netCDF
   C API. The companion program pres_temp_4D_rd.c shows how
   to read the netCDF data file created by this program.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

/* We are writing and reading 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
 * timesteps of data. */
#define NDIMS 4
#define NLAT 6
#define NLON 12
#define LAT_NAME "latitude"
#define LON_NAME "longitude"
#define NREC 2
#define REC_NAME "time"
#define LVL_NAME "level"
#define NLVL 512

/* Names of things. */
#define PRES_NAME "pressure"
#define TEMP_NAME "temperature"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"

/* These are used to construct some example data, and to calculate the values
 * we expect to find. */
#define SAMPLE_PRESSURE 900.0
#define SAMPLE_TEMP 9.0
#define START_LAT 25.0
#define START_LON -125.0

/* For the units attributes. */
#define UNITS "units"
#define PRES_UNITS "hPa"
#define TEMP_UNITS "celsius"
#define LAT_UNITS "degrees_north"
#define LON_UNITS "degrees_east"
#define MAX_ATT_LEN 80

/* This is the name of the data file we will create and read back. */
#define FILE_NAME "pres_temp_4D.nc"


/*----< pres_temp_4D_wr() >--------------------------------------------------*/
int pres_temp_4D_wr(const char *filename,
                    int         format,
                    int         coll_io,
                    MPI_Info    info)
{
    /* IDs for the netCDF file, dimensions, and variables. */
    int nprocs, rank, nerrs = 0, ncid;
    int lon_dimid, lat_dimid, lvl_dimid, rec_dimid;
    int lat_varid, lon_varid, pres_varid, temp_varid;
    int dimids[NDIMS];

    /* The start and count arrays will tell the netCDF library where to
      write our data. */
    MPI_Offset start[NDIMS], count[NDIMS];

    /* Program variables to hold the data we will write out. We will only
      need enough space to hold one timestep of data; one record. */
    float **pres_out; /* [NLVL/nprocs][NLAT][NLON] */
    float **temp_out; /* [NLVL/nprocs][NLAT][NLON] */

    /* These program variables hold the latitudes and longitudes. */
    float lats[NLAT], lons[NLON];

    /* Loop indexes. */
    int lvl, lat, lon, rec, i = 0;

    /* Error handling. */
    int err;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Create some pretend data. If this wasn't an example program, we
     * would have some real data to write, for example, model
     * output.
     */
    for (lat = 0; lat < NLAT; lat++)
        lats[lat] = START_LAT + 5. * lat;
    for (lon = 0; lon < NLON; lon++)
        lons[lon] = START_LON + 5. * lon;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* Create the file. */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
    err = ncmpi_def_dim(ncid, LVL_NAME, NLVL, &lvl_dimid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, LAT_NAME, NLAT, &lat_dimid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, LON_NAME, NLON, &lon_dimid);
    CHECK_ERR
    err = ncmpi_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid);
    CHECK_ERR

    /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
    err = ncmpi_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, &lat_varid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, &lon_varid);
    CHECK_ERR

    /* Assign units attributes to coordinate variables. */
    err = ncmpi_put_att_text(ncid, lat_varid, UNITS,
                             strlen(DEGREES_NORTH), DEGREES_NORTH);
    CHECK_ERR
    err = ncmpi_put_att_text(ncid, lon_varid, UNITS,
                             strlen(DEGREES_EAST), DEGREES_EAST);
    CHECK_ERR

    /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
    dimids[0] = rec_dimid;
    dimids[1] = lvl_dimid;
    dimids[2] = lat_dimid;
    dimids[3] = lon_dimid;

    /* Define the netCDF variables for the pressure and temperature
    * data. */
    err = ncmpi_def_var(ncid, PRES_NAME, NC_FLOAT, NDIMS, dimids, &pres_varid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, TEMP_NAME, NC_FLOAT, NDIMS, dimids, &temp_varid);
    CHECK_ERR

    /* Assign units attributes to the netCDF variables. */
    err = ncmpi_put_att_text(ncid, pres_varid, UNITS,
                             strlen(PRES_UNITS), PRES_UNITS);
    CHECK_ERR
    err = ncmpi_put_att_text(ncid, temp_varid, UNITS,
                             strlen(TEMP_UNITS), TEMP_UNITS);
    CHECK_ERR

    /* End define mode. */
    err = ncmpi_enddef(ncid);
    CHECK_ERR

    err = ncmpi_begin_indep_data(ncid);
    /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
    if (rank == 0) {
        err = ncmpi_put_var_float(ncid, lat_varid, &lats[0]);
        CHECK_ERR
        err = ncmpi_put_var_float(ncid, lon_varid, &lons[0]);
        CHECK_ERR
    }

    if (coll_io) {
        err = ncmpi_end_indep_data(ncid);
        CHECK_ERR
    }

    /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
                    &data[0][0][0]);
     timestep to write.) */
    count[0] = 1;
    count[2] = NLAT;
    count[3] = NLON;
    start[2] = 0;
    start[3] = 0;

    /* divide NLVL dimension among processes */
    count[1] = NLVL / nprocs;
    start[1] = count[1] * rank;
    if (rank < NLVL % nprocs) {
        start[1] += rank;
        count[1]++;
    }
    else {
        start[1] += NLVL % nprocs;
    }
    if (count[1] == 0)
        start[1] = 0;

    /* allocate write buffers */
    pres_out = (float **)malloc(sizeof(float *) * count[1] * 2);
    temp_out = pres_out + count[1];
    if (count[1] > 0) {
        pres_out[0] = (float *)malloc(sizeof(float) * count[1] * 2 * NLAT * NLON);
        temp_out[0] = pres_out[0] + count[1] * NLAT * NLON;
        for (i = 1; i < count[1]; i++) {
            pres_out[i] = pres_out[i - 1] + NLAT * NLON;
            temp_out[i] = temp_out[i - 1] + NLAT * NLON;
        }
    }

    /* initialize write buffers */
    i = (int)start[1] * NLAT * NLON;
    for (lvl = 0; lvl < count[1]; lvl++)
    for (lat = 0; lat < NLAT; lat++)
    for (lon = 0; lon < NLON; lon++) {
        pres_out[lvl][lat * NLON + lon] = SAMPLE_PRESSURE + i;
        temp_out[lvl][lat * NLON + lon] = SAMPLE_TEMP + i++;
    }

    /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */

    for (rec = 0; rec < NREC; rec++) {
        start[0] = rec;
        if (coll_io)
            err = ncmpi_put_vara_float_all(ncid, pres_varid, start, count, &pres_out[0][0]);
        else
            err = ncmpi_put_vara_float(ncid, pres_varid, start, count, &pres_out[0][0]);
        CHECK_ERR
        if (coll_io)
            err = ncmpi_put_vara_float_all(ncid, temp_varid, start, count, &temp_out[0][0]);
        else
            err = ncmpi_put_vara_float(ncid, temp_varid, start, count, &temp_out[0][0]);
        CHECK_ERR
    }

    /* Close the file. */
    err = ncmpi_close(ncid);
    CHECK_ERR

    if (count[1] > 0)
        free(pres_out[0]);
    free(pres_out);

    return (nerrs > 0);
}

/*----< pres_temp_4D_rd() >--------------------------------------------------*/
int pres_temp_4D_rd(const char *filename,
                    int         coll_io,
                    MPI_Info    info)
{
    int rank, nprocs, ncid, pres_varid, temp_varid, lat_varid, lon_varid;

    /* The start and count arrays will tell the netCDF library where to
      read our data. */
    MPI_Offset start[NDIMS], count[NDIMS];

    /* Program variables to hold the data we will read. We will only
      need enough space to hold one timestep of data; one record. */
    float **pres_in = NULL; /* [NLVL/nprocs][NLAT][NLON] */
    float **temp_in = NULL; /* [NLVL/nprocs][NLAT][NLON] */

    /* These program variables hold the latitudes and longitudes. */
    float lats[NLAT], lons[NLON];

    /* Loop indexes. */
    int lvl, lat, lon, rec, i = 0;

    /* Error handling. */
    int err, nerrs = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Open the file. */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    if (err != NC_NOERR) { /* fatal error */
        if (rank == 0)
            fprintf(stderr,"Error: failed to open file %s (%s)\n",filename,ncmpi_strerror(err));
        goto err_out;
    }

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* Get the varids of the latitude and longitude coordinate variables. */
    err = ncmpi_inq_varid(ncid, LAT_NAME, &lat_varid);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, LON_NAME, &lon_varid);
    CHECK_ERR

    /* Read the coordinate variable data. */
    memset(lats, 0, sizeof(float) * NLAT);
    memset(lons, 0, sizeof(float) * NLON);
    if (coll_io) {
        err = ncmpi_get_var_float_all(ncid, lat_varid, &lats[0]);
        CHECK_ERR
        err = ncmpi_get_var_float_all(ncid, lon_varid, &lons[0]);
        CHECK_ERR
    }
    else {
        err = ncmpi_get_var_float(ncid, lat_varid, &lats[0]);
        CHECK_ERR
        err = ncmpi_get_var_float(ncid, lon_varid, &lons[0]);
        CHECK_ERR
    }

    /* Check the coordinate variable data. */
    for (lat = 0; lat < NLAT; lat++) {
        float exp =  START_LAT + 5. * lat;
        if (lats[lat] != exp) {
            printf("\nError at line %d in %s: %s[%d] expect %.1f but got %.1f\n",
                   __LINE__, __FILE__, LAT_NAME, lat, exp, lats[lat]);
            nerrs++;
            break;
        }
    }
    for (lon = 0; lon < NLON; lon++) {
        float exp =  START_LON + 5. * lon;
        if (lons[lon] != START_LON + 5. * lon) {
            printf("\nError at line %d in %s: %s[%d] expect %.1f but got %.1f\n",
                   __LINE__, __FILE__, LON_NAME, lon, exp, lons[lon]);
            nerrs++;
            break;
        }
    }

    /* Get the varids of the pressure and temperature netCDF variables. */
    err = ncmpi_inq_varid(ncid, PRES_NAME, &pres_varid);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, TEMP_NAME, &temp_varid);
    CHECK_ERR

    /* Read the data. Since we know the contents of the file we know that the
     * data arrays in this program are the correct size to hold one timestep.
     */
    count[0] = 1;
    count[2] = NLAT;
    count[3] = NLON;
    start[2] = 0;
    start[3] = 0;

    /* divide NLVL dimension among processes */
    count[1] = NLVL / nprocs;
    start[1] = count[1] * rank;
    if (rank < NLVL % nprocs) {
        start[1] += rank;
        count[1]++;
    }
    else {
        start[1] += NLVL % nprocs;
    }
    if (count[1] == 0)
        start[1] = 0;

    /* allocate read buffers */
    pres_in = (float **)malloc(sizeof(float *) * count[1] * 2);
    temp_in = pres_in + count[1];
    if (count[1] > 0) {
        pres_in[0] = (float *)malloc(sizeof(float) * count[1] * 2 * NLAT * NLON);
        temp_in[0] = pres_in[0] + count[1] * NLAT * NLON;
        for (i = 1; i < count[1]; i++) {
            pres_in[i] = pres_in[i - 1] + NLAT * NLON;
            temp_in[i] = temp_in[i - 1] + NLAT * NLON;
        }
    }

    /* Read and check one record at a time. */
    for (rec = 0; rec < NREC; rec++) {
        start[0] = rec;
        if (coll_io)
            err = ncmpi_get_vara_float_all(ncid, pres_varid, start, count, &pres_in[0][0]);
        else
            err = ncmpi_get_vara_float(ncid, pres_varid, start, count, &pres_in[0][0]);
        CHECK_ERR
        if (coll_io)
            err = ncmpi_get_vara_float_all(ncid, temp_varid, start, count, &temp_in[0][0]);
        else
            err = ncmpi_get_vara_float(ncid, temp_varid, start, count, &temp_in[0][0]);
        CHECK_ERR

        /* Check the data. */
        i = (int)start[1] * NLAT * NLON;
        for (lvl = 0; lvl < count[1]; lvl++)
        for (lat = 0; lat < NLAT; lat++)
        for (lon = 0; lon < NLON; lon++) {
            float exp = SAMPLE_PRESSURE + i;
            int indx = lat * NLON + lon;
            if (pres_in[lvl][indx] != exp) {
                printf("\nError at line %d in %s: %s[%d][%d][%d][%d] expect %.1f but got %.1f\n",
                __LINE__, __FILE__, PRES_NAME, rec, lvl, lat, lon, exp, pres_in[lvl][indx]);
                nerrs++;
                goto fn_exit;
            }
            exp = SAMPLE_TEMP + i;
            if (temp_in[lvl][indx] != exp) {
                printf("\nError at line %d in %s: %s[%d][%d][%d][%d] expect %.1f but got %.1f\n",
                __LINE__, __FILE__, TEMP_NAME, rec, lvl, lat, lon, exp, temp_in[lvl][indx]);
                nerrs++;
                goto fn_exit;
            }
            i++;
        }

fn_exit:
        MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (nerrs > 0) break;
    } /* next record */

    /* Close the file. */
    err = ncmpi_close(ncid);
    CHECK_ERR

    if (pres_in != NULL) {
        if (pres_in[0] != NULL)
            free(pres_in[0]);
        free(pres_in);
    }

err_out:
    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int nerrs;

    nerrs = pres_temp_4D_wr(out_path, format, coll_io, info);
    if (nerrs > 0) return nerrs;

    MPI_Barrier(MPI_COMM_WORLD);

    nerrs = pres_temp_4D_rd(out_path, coll_io, info);
    if (nerrs > 0) return nerrs;

    return 0;
}

#if 0
static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-q] quiet mode\n"
    "       [-k] Keep output files (default: no)\n"
    "       [filename]: output netCDF file name (default: %s)\n";
    fprintf(stderr, help, basename(argv0), FILE_NAME);
}

int tst_main(int    argc,
             char **argv,
             char  *msg,
             int   (*tst_body)(const char*, int, int, MPI_Info))
{
    extern int optind;
    extern char *optarg;
    char filename[256], *ptr;

    /* IDs for the netCDF file, dimensions, and variables. */
    int nprocs, rank, err, nerrs=0, keep_files, quiet, coll_io;
    int i, a, d, r, m, b;
#ifdef ENABLE_NETCDF4
    int num_fmts = 4;
    int format[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_NETCDF4_CLASSIC, NC_FORMAT_64BIT_DATA};
#else
    int num_fmts = 3;
    int format[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
#endif

    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    keep_files = 0;
    quiet = 0;

    while ((i = getopt(argc, argv, "hqk")) != EOF)
        switch(i) {
            case 'q':
                quiet = 1;
                break;
            case 'k':
                keep_files = 1;
                break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for %s", basename(argv[0]), msg);
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    char cmd_opts[64];
    sprintf(cmd_opts, "Rank %d: ncmpidiff", rank);
    ptr = strrchr(filename, '.');
    if (ptr != NULL) *ptr = '\0';

    MPI_Info_create(&info);

    for (i=0; i<num_fmts; i++) {
        char path[512], ext[16], *base_file;

        base_file = NULL;
        sprintf(ext, "nc%d", format[i]);

        for (a=0; a<2; a++) {
        for (d=0; d<2; d++) {
        for (r=0; r<2; r++) {
        for (m=0; m<2; m++) {
        for (b=0; b<2; b++) {

            sprintf(path, "%s.%s", filename, ext);

            if (a == 0) {
                MPI_Info_set(info, "nc_num_aggrs_per_node", "0");
                strcat(path, ".noina");
            } else {
                MPI_Info_set(info, "nc_num_aggrs_per_node", "2");
                strcat(path, ".ina");
            }

            if (d == 0) {
                MPI_Info_set(info, "nc_pncio", "enable");
                strcat(path, ".pncio");
            } else {
                MPI_Info_set(info, "nc_pncio", "disable");
                strcat(path, ".mpio");
            }

            if (r == 0) {
                MPI_Info_set(info, "romio_no_indep_rw", "false");
                strcat(path, ".noindep");
            } else {
                MPI_Info_set(info, "romio_no_indep_rw", "true");
                strcat(path, ".indep");
            }

            if (m == 0) {
                MPI_Info_set(info, "nc_data_move_chunk_size", "100");
                strcat(path, ".chunk100");
            } else {
                MPI_Info_set(info, "nc_data_move_chunk_size", "1048576");
                strcat(path, ".chunk1M");
            }

            if (b == 0) {
                MPI_Info_set(info, "nc_burst_buf", "disable");
                strcat(path, ".nobb");
            }
            else {
#ifdef ENABLE_BURST_BUFFER
                MPI_Info_set(info, "nc_burst_buf", "enable");
                MPI_Info_set(info, "nc_burst_buf_dirname", TESTOUTDIR);
                MPI_Info_set(info, "nc_burst_buf_overwrite", "enable");
                strcat(path, ".bb");
#else
                continue;
#endif
            }

            for (coll_io=0; coll_io<2; coll_io++) {

                /* NetCDF4 does not allow to extend number of record numbers in
                 * independent data mode. NC_ECANTEXTEND will be returned.
                 */
                if (format[i] == NC_FORMAT_NETCDF4_CLASSIC && coll_io == 0)
                    continue;

                nerrs = tst_body(path, format[i], coll_io, info);
                if (nerrs != NC_NOERR) goto err_out;
            }

            /* run ncmpidiff to compare output files */
            if (base_file == NULL) /* skip first file */
                base_file = strdup(path);
            else if (strcmp(base_file, path)) {
                int check_header=1, check_entire_file=1, first_diff=1;

                /* ncmpidiff does nott support netCDF4 files */
                if (format[i] == NC_FORMAT_NETCDF4_CLASSIC) {
                    if (!keep_files) ncmpi_delete(path, MPI_INFO_NULL);
                    continue;
                }

                if (!quiet && rank == 0)
                    printf("ncmpidiff %-60s %s a=%d d=%d r=%d m=%d b=%d\n",
                           path, base_file,a,d,r,m,b);

#ifdef MIMIC_LUSTRE
                /* use a larger stripe size when running ncmpidiff */
                setenv("MIMIC_STRIPE_SIZE", "1048576", 1);
#else
                unsetenv("MIMIC_STRIPE_SIZE");
#endif
                /* running ncmpidiff also validates the file header */
                nerrs = ncmpidiff_core(path, base_file, MPI_COMM_WORLD,
                                       MPI_INFO_NULL, 0, quiet, check_header,
                                       0, check_entire_file, 0, NULL, 0,
                                       first_diff, cmd_opts, 0, 0);
                if (nerrs != 0) goto err_out;
                if (!keep_files) ncmpi_delete(path, MPI_INFO_NULL);
            }
        } /* loop b */
        } /* loop m */
        } /* loop r */
        } /* loop d */
        } /* loop a */

        if (base_file != NULL) {
            if (!keep_files) ncmpi_delete(base_file, MPI_INFO_NULL);
            free(base_file);
        }
    }
    MPI_Info_free(&info);

    /* check if there is any malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

err_out:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs)
            printf(FAIL_STR, nerrs);
        else
            printf(PASS_STR);
    }

    MPI_Finalize();

    return (nerrs > 0);
}
#endif

int main(int argc, char **argv) {

#ifdef ENABLE_NETCDF4
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_NETCDF4_CLASSIC, NC_FORMAT_64BIT_DATA};
#else
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
#endif

    loop_opts opt;

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1;
    opt.drv      = 1;
    opt.ind      = 1;
    opt.chk      = 1;
    opt.bb       = 1;
    opt.mod      = 1;
    opt.hdr_diff = 1;
    opt.var_diff = 1;

    return tst_main(argc, argv, "write/read netCDF file", opt, test_io);
}

