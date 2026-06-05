/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
   This is an example which reads some 4D pressure and
   temperatures. The data file read by this program is produced by the
   companion program pres_temp_4D_wr.c. It is intended to illustrate
   the use of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id$
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>
#include <mpi.h>
#include <testutils.h>

#include "pres_temp_4D.h"

int main(int argc, char **argv)
{
    int rank, nprocs, ncid, pres_varid, temp_varid;
    int lat_varid, lon_varid;

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

    char *filename = FILE_NAME;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 3) {
        if (!rank)
            printf("Usage: %s [filename]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    if (argc > 1) filename = argv[1];

    /* Open the file. */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) { /* fatal error */
        if (rank == 0)
            fprintf(stderr,"Error: failed to open file %s (%s)\n",filename,ncmpi_strerror(err));
        goto err_out;
    }

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
        int format;
        err = ncmpi_inq_format(ncid, &format); CHECK_ERR
        if (format == NC_FORMAT_NETCDF4)
            sprintf(cmd_str, "*** TESTING C   %s for reading NetCDF-4 file", basename(argv[0]));
        else
            sprintf(cmd_str, "*** TESTING C   %s for reading classic file", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    /* Get the varids of the latitude and longitude coordinate variables. */
    err = ncmpi_inq_varid(ncid, LAT_NAME, &lat_varid);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, LON_NAME, &lon_varid);
    CHECK_ERR

    /* Read the coordinate variable data. */
    memset(lats, 0, sizeof(float) * NLAT);
    memset(lons, 0, sizeof(float) * NLON);
    err = ncmpi_get_var_float_all(ncid, lat_varid, &lats[0]);
    CHECK_ERR
    err = ncmpi_get_var_float_all(ncid, lon_varid, &lons[0]);
    CHECK_ERR

    /* Check the coordinate variable data. */
    for (lat = 0; lat < NLAT; lat++) {
        float exp =  START_LAT + 5. * lat;
        if (lats[lat] != exp) {
            printf("\nError at line %d in %s: lats[%d] expect %.1f but got %.1f\n",
            __LINE__, __FILE__, lat, exp, lats[lat]);
            nerrs++;
            break;
        }
    }
    for (lon = 0; lon < NLON; lon++) {
        float exp =  START_LON + 5. * lon;
        if (lons[lon] != START_LON + 5. * lon) {
            printf("\nError at line %d in %s: lons[%d] expect %.1f but got %.1f\n",
            __LINE__, __FILE__, lon, exp, lons[lon]);
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
        err = ncmpi_get_vara_float_all(ncid, pres_varid, start, count, &pres_in[0][0]);
        CHECK_ERR
        err = ncmpi_get_vara_float_all(ncid, temp_varid, start, count, &temp_in[0][0]);
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
