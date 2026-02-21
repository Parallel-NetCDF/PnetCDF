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
                goto loop_exit;
            }
            exp = SAMPLE_TEMP + i;
            if (temp_in[lvl][indx] != exp) {
                printf("\nError at line %d in %s: %s[%d][%d][%d][%d] expect %.1f but got %.1f\n",
                __LINE__, __FILE__, TEMP_NAME, rec, lvl, lat, lon, exp, temp_in[lvl][indx]);
                nerrs++;
                goto loop_exit;
            }
            i++;
        }

loop_exit:
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

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "write/read netCDF file", opt, test_io);

    MPI_Finalize();

    return err;
}

