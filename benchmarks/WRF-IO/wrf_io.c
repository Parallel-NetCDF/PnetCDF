/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program evaluates the file write and read performance of WRF (Weather
 * Research and Forecast Model, https://github.com/wrf-model/WRF) developed at
 * NCAR. It's data partitioning pattern is a 2D block-block checkerboard
 * pattern, along the longitude and latitude. This benchmark program reads a
 * CDL header file in text format (i.e. a netCDF file obtained from "ncmpidump
 * -h") and creates a new file with the same metadata.
 *
 * The provided CDL file "wrf_header.txt" contains 202 variables. Among them,
 * 30 are 1D variables. 25 are 2D variables, 122 are 3D, and 25 are 4D. Only 3D
 * and 4D variables are partitioned. 1D and 2D variables are only written by
 * process of rank 0.
 *
 * The grid size can be set by command-line options of -l and -w.
 * -l to set size for dimension "south_north"
 * -w to set size for dimension "west_east"
 *
 * To compile:
 *    % mpicc -O2 wrf_io.c -o wrf_io -lpnetcdf
 *
 * An example of run command:
 *    % mpiexec -n 4 ./wrf_io -y 100 -x 100 -n 2 -i wrf_header.txt -w ./wrf_io.nc
 *
 *    -----------------------------------------------------------
 *    ---- WRF-IO write benchmark ----
 *    Output NetCDF file name:        ./wrf_io.nc
 *    Number of MPI processes:        4
 *    MPI processes arranged to 2D:   2 x 2
 *    Grid size longitude x latitude: 100 x 100
 *    Total number of variables:      202
 *    Number of time records:         2
 *    Total write amount:             68696794 B
 *                                    65.51 MiB
 *                                    0.06 GiB
 *    Max open-to-close time:         0.0978 sec
 *    Max define metadata time:       0.0085 sec
 *    Max bput posting  time:         0.0148 sec
 *    Max wait_all      time:         0.0718 sec
 *    Write bandwidth:                670.08 MiB/s
 *                                    0.65 GiB/s
 *    -----------------------------------------------------------
 *    MPI-IO hint striping_factor:        0
 *    MPI-IO hint striping_unit:          0
 *    MPI-IO hint cb_buffer_size:         16777216
 *    MPI-IO hint aggr_list:              16777216
 *    MPI-IO hint cb_nodes:               1
 *    PnetCDF hint nc_num_aggrs_per_node: 0
 *    -----------------------------------------------------------
 *    Max heap memory allocated by PnetCDF internally is 8.54 MiB
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy(), strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

static int verbose, debug;

#define CHECK_ERR(name) {                               \
    if (err != NC_NOERR) {                              \
        printf("Error at line=%d: name=%s error=%s\n",  \
               __LINE__, name, ncmpi_strerror(err));    \
        goto err_out;                                   \
    }                                                   \
}

#define HINT ((flag)?(value):("NOT SET"))

typedef struct {
    int varid;
    char *name;
    int ndims;
    int xtype;
    int *dimids;
    MPI_Offset nelems;
    MPI_Offset start[4];
    MPI_Offset count[4];
    void *buf;
} WRF_VAR;

static
int construct_vars(int         hid, /* CDL header ID */
                   WRF_VAR    *vars,
                   MPI_Offset *longitude,
                   MPI_Offset *latitude,
                   int         psizes[2],
                   MPI_Offset *buf_size)
{
    char *name;
    int i, err=NC_NOERR, nprocs, rank;
    int ndims, *dimids, nvars, my_rank_y, my_rank_x;
    MPI_Offset dimlen, my_start_y, my_start_x, my_count_y, my_count_x;
    nc_type xtype;

    /* determine whether to use longitude, latitude set at command line */

    /* retrieve the number of dimensions defined in the CDL file */
    err = cdl_hdr_inq_ndims(hid, &ndims);
    CHECK_ERR("cdl_hdr_inq_ndims")

    if (*longitude <= 0) {
        for (i=0; i<ndims; i++) {
            /* retrieve metadata of dimension i */
            err = cdl_hdr_inq_dim(hid, i, &name, &dimlen);
            CHECK_ERR("cdl_hdr_inq_dim")

            /* Note south_north must be defined before south_north_stag */
            if (!strcmp(name, "south_north")) {
                *longitude = dimlen;
                break;
            }
        }
    }
    if (*latitude <= 0) {
        for (i=0; i<ndims; i++) {
            /* retrieve metadata of dimension i */
            err = cdl_hdr_inq_dim(hid, i, &name, &dimlen);
            CHECK_ERR("cdl_hdr_inq_dim")

            /* Note west_east must be defined before west_east_stag */
            if (!strcmp(name, "west_east")) {
                *latitude = dimlen;
                break;
            }
        }
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    my_rank_y = rank / psizes[1];
    my_rank_x = rank % psizes[1];

    my_count_y = *longitude / psizes[0];
    my_start_y = my_count_y * my_rank_y;
    if (my_rank_y < *longitude % psizes[0]) {
        my_start_y += my_rank_y;
        my_count_y++;
    }
    else {
        my_start_y += *longitude % psizes[0];
    }

    my_count_x = *latitude / psizes[1];
    my_start_x = my_count_x * my_rank_x;
    if (my_rank_x < *latitude % psizes[1]) {
        my_start_x += my_rank_x;
        my_count_x++;
    }
    else {
        my_start_x += *latitude % psizes[1];
    }

    if (debug) {
        printf("%2d: rank (%2d, %2d) start %4lld %4lld count %4lld %4lld\n",
               rank, my_rank_y, my_rank_x, my_start_y, my_start_x, my_count_y, my_count_x);
        fflush(stdout);
    }

    *buf_size = 0;

    /* retrieve number of variables defined in the CDL file */
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")
    if (debug) printf("var: nvars %d\n", nvars);

    for (i=0; i<nvars; i++, vars++) {
        char *dname1, *dname2, *dname3;
        MPI_Offset dim1, dim2, dim3;

        /* retrieve metadata of variable i defined in the CDL file */
        err = cdl_hdr_inq_var(hid, i, &name, &xtype, &ndims, &dimids);
        CHECK_ERR(name)

        vars->name = name;
        vars->ndims = ndims;
        vars->xtype = xtype;
        vars->dimids = dimids;
        vars->start[0] = 0;    /* time dimension */
        vars->count[0] = 1;    /* time dimension */

        if (ndims == 1) {
            vars->nelems = (rank == 0) ? 1 : 0;
        }
        else if (ndims == 2) {
            err = cdl_hdr_inq_dim(hid, dimids[1], NULL, &dim1);
            CHECK_ERR("cdl_hdr_inq_dim")
            vars->nelems = (rank == 0) ? dim1 : 0;
            vars->count[1] = dim1; /* dimension dim1 is not partitioned */
            vars->start[1] = 0;    /* dimension dim1 is not partitioned */
        }
        else if (ndims == 3) {
            err = cdl_hdr_inq_dim(hid, dimids[1], &dname1, &dim1);
            CHECK_ERR("cdl_hdr_inq_dim")
            err = cdl_hdr_inq_dim(hid, dimids[2], &dname2, &dim2);
            CHECK_ERR("cdl_hdr_inq_dim")
            vars->start[1] = my_start_y;
            vars->count[1] = my_count_y;
            vars->start[2] = my_start_x;
            vars->count[2] = my_count_x;
            if (!strcmp(dname1, "south_north_stag") && my_rank_y == psizes[0]-1)
                vars->count[1]++;
            if (!strcmp(dname2, "west_east_stag") && my_rank_x == psizes[1]-1)
                vars->count[2]++;
            vars->nelems = vars->count[1] *  vars->count[2];
        }
        else if (ndims == 4) {
            err = cdl_hdr_inq_dim(hid, dimids[1], NULL, &dim1);
            CHECK_ERR("cdl_hdr_inq_dim")
            err = cdl_hdr_inq_dim(hid, dimids[2], &dname2, &dim2);
            CHECK_ERR("cdl_hdr_inq_dim")
            err = cdl_hdr_inq_dim(hid, dimids[3], &dname3, &dim3);
            CHECK_ERR("cdl_hdr_inq_dim")
            vars->start[1] = 0;    /* this dimension is not partitioned */
            vars->count[1] = dim1; /* this dimension is not partitioned */
            vars->start[2] = my_start_y;
            vars->count[2] = my_count_y;
            vars->start[3] = my_start_x;
            vars->count[3] = my_count_x;
            if (!strcmp(dname2, "south_north_stag") && my_rank_y == psizes[0]-1)
                vars->count[2]++;
            if (!strcmp(dname3, "west_east_stag") && my_rank_x == psizes[1]-1)
                vars->count[3]++;
            vars->nelems = dim1 * vars->count[2] *  vars->count[3];
        }

        if (xtype == NC_FLOAT)
            *buf_size += sizeof(float) * vars->nelems;
        else if (xtype == NC_INT)
            *buf_size += sizeof(int) * vars->nelems;
        else if (xtype == NC_CHAR)
            *buf_size += vars->nelems;
    }

err_out:
    return err;
}

static
int inquire_vars(int         ncid,
                 WRF_VAR    *vars,
                 int         psizes[2],
                 MPI_Offset  longitude,
                 MPI_Offset  latitude,
                 MPI_Offset *buf_size)
{
    int i, j, err=NC_NOERR, nprocs, rank;
    int ndims, dimids[4], nvars, my_rank_y, my_rank_x;
    MPI_Offset my_start_y, my_start_x, my_count_y, my_count_x;
    nc_type xtype;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    my_rank_y = rank / psizes[1];
    my_rank_x = rank % psizes[1];

    my_count_y = longitude / psizes[0];
    my_start_y = my_count_y * my_rank_y;
    if (my_rank_y < longitude % psizes[0]) {
        my_start_y += my_rank_y;
        my_count_y++;
    }
    else {
        my_start_y += longitude % psizes[0];
    }

    my_count_x = latitude / psizes[1];
    my_start_x = my_count_x * my_rank_x;
    if (my_rank_x < latitude % psizes[1]) {
        my_start_x += my_rank_x;
        my_count_x++;
    }
    else {
        my_start_x += latitude % psizes[1];
    }

    if (debug) {
        printf("%2d: rank (%2d, %2d) start %4lld %4lld count %4lld %4lld\n",
               rank, my_rank_y, my_rank_x, my_start_y, my_start_x, my_count_y, my_count_x);
        fflush(stdout);
    }

    *buf_size = 0;

    /* inquire the number of fix-sized and record variables */
    err = ncmpi_inq_nvars(ncid, &nvars);
    CHECK_ERR("ncmpi_inq_nvars")

    for (i=0; i<nvars; i++, vars++) {
        char name[64], dname0[64], dname1[64], dname2[64], dname3[64];
        MPI_Offset dim0, dim1, dim2, dim3;
        int nattrs;

        err = ncmpi_inq_var(ncid, i, name, &xtype, &ndims, dimids, &nattrs);
        CHECK_ERR(name)

        if (ndims > 0) {
            err = ncmpi_inq_dimlen(ncid, dimids[0], &dim0);
            CHECK_ERR("ncmpi_inq_dimlen")
            err = ncmpi_inq_dimname(ncid, dimids[0], dname0);
            CHECK_ERR("ncmpi_inq_dimname")
        }
        if (ndims > 1) {
            err = ncmpi_inq_dimlen(ncid, dimids[1], &dim1);
            CHECK_ERR("ncmpi_inq_dimlen")
            err = ncmpi_inq_dimname(ncid, dimids[1], dname1);
            CHECK_ERR("ncmpi_inq_dimname")
        }
        if (ndims > 2) {
            err = ncmpi_inq_dimlen(ncid, dimids[2], &dim2);
            CHECK_ERR("ncmpi_inq_dimlen")
            err = ncmpi_inq_dimname(ncid, dimids[2], dname2);
            CHECK_ERR("ncmpi_inq_dimname")
        }
        if (ndims > 3) {
            err = ncmpi_inq_dimlen(ncid, dimids[3], &dim3);
            CHECK_ERR("ncmpi_inq_dimlen")
            err = ncmpi_inq_dimname(ncid, dimids[3], dname3);
            CHECK_ERR("ncmpi_inq_dimname")
        }

        vars->varid  = i;
        vars->name   = strdup(name);
        vars->ndims  = ndims;
        vars->xtype  = xtype;
        vars->nelems = 0;
        vars->dimids = (int*) malloc(sizeof(int) * ndims);
        for (j=0; j<ndims; j++) vars->dimids[j] = dimids[j];

        vars->start[0] = 0;    /* time dimension */
        vars->count[0] = 1;    /* time dimension */

        /* In WRF, the first dimension is always NC_UNLIMITED */
        if (ndims == 1)
            vars->nelems = (rank == 0) ? 1 : 0;
        else if (ndims == 2) {
            vars->nelems = (rank == 0) ? dim1 : 0;
            vars->count[1] = dim1; /* dimension dim1 is not partitioned */
            vars->start[1] = 0;    /* dimension dim1 is not partitioned */
        }
        else if (ndims == 3) {
            vars->start[1] = my_start_y;
            vars->count[1] = my_count_y;
            vars->start[2] = my_start_x;
            vars->count[2] = my_count_x;
            if (!strcmp(dname1, "south_north_stag") && my_rank_y == psizes[0]-1)
                vars->count[1]++;
            if (!strcmp(dname2, "west_east_stag") && my_rank_x == psizes[1]-1)
                vars->count[2]++;
            vars->nelems = vars->count[1] *  vars->count[2];
        }
        else if (ndims == 4) {
            vars->start[1] = 0;    /* this dimension is not partitioned */
            vars->count[1] = dim1; /* this dimension is not partitioned */
            vars->start[2] = my_start_y;
            vars->count[2] = my_count_y;
            vars->start[3] = my_start_x;
            vars->count[3] = my_count_x;
            if (!strcmp(dname2, "south_north_stag") && my_rank_y == psizes[0]-1)
                vars->count[2]++;
            if (!strcmp(dname3, "west_east_stag") && my_rank_x == psizes[1]-1)
                vars->count[3]++;
            vars->nelems = dim1 * vars->count[2] *  vars->count[3];
        }

        if (xtype == NC_FLOAT)
            *buf_size += sizeof(float) * vars->nelems;
        else if (xtype == NC_INT)
            *buf_size += sizeof(int) * vars->nelems;
        else if (xtype == NC_CHAR)
            *buf_size += vars->nelems;
    }

err_out:
    return err;
}

static
int def_dims_vars(int         ncid,
                  int         hid, /* CDL header ID */
                  MPI_Offset *longitude,
                  MPI_Offset *latitude,
                  WRF_VAR    *vars)
{
    char *name;
    void *value;
    int i, j, err=NC_NOERR, ndims, nvars, nattrs;
    MPI_Offset size, nelems;
    nc_type xtype;

    /* define dimensions */

    /* retrieve the number of dimensions defined in the CDL file */
    err = cdl_hdr_inq_ndims(hid, &ndims);
    CHECK_ERR("cdl_hdr_inq_ndims")
    if (debug) printf("dim: ndims %d\n", ndims);

    for (i=0; i<ndims; i++) {
        int dimid;

        /* retrieve metadata of dimension i */
        err = cdl_hdr_inq_dim(hid, i, &name, &size);
        CHECK_ERR("cdl_hdr_inq_dim")
        if (debug) printf("\t name %s size %lld\n",name, size);

             if (!strcmp(name, "south_north"))      size = *longitude;
        else if (!strcmp(name, "south_north_stag")) size = *longitude + 1;
        else if (!strcmp(name, "west_east"))        size = *latitude;
        else if (!strcmp(name, "west_east_stag"))   size = *latitude + 1;

        err = ncmpi_def_dim(ncid, name, size, &dimid);
        CHECK_ERR("ncmpi_def_dim")
    }

    /* retrieve number of variables defined in the CDL file */
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")

    /* define variables */
    for (i=0; i<nvars; i++, vars++) {
        err = ncmpi_def_var(ncid, vars->name, vars->xtype, vars->ndims,
                            vars->dimids, &vars->varid);
        CHECK_ERR(name)

        /* define local attributes */

        /* retrieve metadata of attribute j associated with variable i */
        err = cdl_hdr_inq_nattrs(hid, i, &nattrs);
        CHECK_ERR("cdl_hdr_inq_nattrs")

        if (debug) {
            printf("\t name %s type %d ndims %d nattr %d\n",
                          name, vars->xtype, ndims, nattrs);
            for (j=0; j<ndims; j++)
                printf("\t\tdimid %d\n",vars->dimids[j]);
        }

        for (j=0; j<nattrs; j++) {
            /* retrieve metadata of attribute j associated with variable i */
            err = cdl_hdr_inq_attr(hid, i, j, &name, &xtype, &nelems, &value);
            CHECK_ERR("cdl_hdr_inq_attr")
            if (debug) {
                if (xtype == NC_CHAR)
                    printf("\t\tattr %s type %d nelems %lld (%s)\n",
                            name, xtype,nelems,(char*)value);
                else
                    printf("\t\tattr %s type %d nelems %lld\n",
                           name, xtype, nelems);
            }

            err = ncmpi_put_att(ncid, vars->varid, name, xtype, nelems, value);
            CHECK_ERR("ncmpi_put_att")
        }
    }

    /* define global attributes */

    /* retrieve the number of global attributes */
    err = cdl_hdr_inq_nattrs(hid, NC_GLOBAL, &nattrs);
    CHECK_ERR("cdl_hdr_inq_nattrs")
    if (debug) printf("global attrs: nattrs %d\n", nattrs);

    for (i=0; i<nattrs; i++) {
        /* retrieve metadata of global attribute i */
        err = cdl_hdr_inq_attr(hid, NC_GLOBAL, i, &name, &xtype, &nelems, &value);
        CHECK_ERR("cdl_hdr_inq_attr")
        if (debug) {
            if (xtype == NC_CHAR)
                printf("\t name %s type %d nelems %lld (%s)\n",
                        name, xtype, nelems,(char*)value);
            else
                printf("\t name %s type %d nelems %lld\n",
                        name, xtype, nelems);
        }

        err = ncmpi_put_att(ncid, NC_GLOBAL, name, xtype, nelems, value);
        CHECK_ERR("ncmpi_put_att")
    }

err_out:
    return err;
}

static
int wrf_w_benchmark(char       *out_file,
                    int         hid,  /* CDL header ID */
                    int         psizes[2],
                    MPI_Offset  longitude,
                    MPI_Offset  latitude,
                    int         ntimes,
                    MPI_Info    info)
{
    int i, j, err=NC_NOERR, nprocs, rank;
    int cmode, ncid, nvars;
    double timing[4], max_t[4];
    MPI_Offset buf_size, w_size, sum_w_size;
    MPI_Info info_used;
    WRF_VAR *vars=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* allocate space for storing variable metadata */

    /* retrieve number of variables defined in the CDL file */
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")

    vars = (WRF_VAR*) malloc(sizeof(WRF_VAR) * nvars);

    /* populate variable metadata */
    err = construct_vars(hid, vars, &longitude, &latitude, psizes, &buf_size);
    CHECK_ERR("construct_vars")

    if (debug) {
        printf("%2d: buf_size %lld\n", rank, buf_size);
        fflush(stdout);
    }

    /* allocate and initialize write buffers */
    MPI_Offset mem_alloc;
    if (debug) mem_alloc = 0;

    for (i=0; i<nvars; i++) {
        if (vars[i].nelems == 0) continue;

        if (vars[i].xtype == NC_FLOAT) {
            float *flt_ptr = (float*) malloc(sizeof(float) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                flt_ptr[j] = rank + j;
            vars[i].buf = (void*) flt_ptr;
            if (debug) mem_alloc += sizeof(float) * vars[i].nelems;
        }
        else if (vars[i].xtype == NC_INT) {
            int *int_ptr = (int*) malloc(sizeof(int) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                int_ptr[j] = rank + j;
            vars[i].buf = (void*) int_ptr;
            if (debug) mem_alloc += sizeof(int) * vars[i].nelems;
        }
        else if (vars[i].xtype == NC_CHAR) {
            char *str_ptr = (char*) malloc(vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                str_ptr[j] = '0' + rank + j;
            vars[i].buf = (void*) str_ptr;
            if (debug) mem_alloc += vars[i].nelems;
        }
    }

    if (debug) {
        if (mem_alloc != buf_size)
            printf("%2d: mem_alloc %lld != buf_size %lld\n", rank, mem_alloc,
                   buf_size);
        fflush(stdout);
        assert(mem_alloc == buf_size);
    }

    /* start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();

    /* create the output file */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, out_file, cmode, info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line=%d: creating file %s (%s)\n",
               __LINE__, out_file, ncmpi_strerror(err));
        goto err_out;
    }

    /* define dimension, variables, and attributes */
    err = def_dims_vars(ncid, hid, &longitude, &latitude, vars);
    CHECK_ERR("def_dims_vars")

    /* exit metadata define mode */
    err = ncmpi_enddef(ncid);
    CHECK_ERR("ncmpi_enddef")

    /* tell PnetCDF how much space for storing pending write requests before
     * flushing it.
     */
    err = ncmpi_buffer_attach(ncid, buf_size);
    CHECK_ERR("ncmpi_buffer_attach")

    timing[1] = MPI_Wtime() - timing[0];
    timing[2] = timing[3] = 0;

    /* ntimes is the number of records to be written */
    for (j=0; j<ntimes; j++) {
        double start_t, end_t;
        start_t = MPI_Wtime();

        if (debug && rank == 0) {
            printf("Writing record %d\n",j);
            fflush(stdout);
        }

        for (i=0; i<nvars; i++) {
            if (vars[i].nelems == 0) continue;

            /* set record ID */
            vars[i].start[0] = j;

            if (vars[i].xtype == NC_FLOAT)
                err = ncmpi_bput_vara_float(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].xtype == NC_INT)
                err = ncmpi_bput_vara_int(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].xtype == NC_CHAR)
                err = ncmpi_bput_vara_text(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            CHECK_ERR(vars[i].name)

        }
        end_t = MPI_Wtime();
        timing[2] += end_t - start_t;
        start_t = end_t;

        if (debug && rank == 0) {
            printf("Flush write requests at end of iteration j=%d\n",j);
            fflush(stdout);
        }

        /* flush all nonblocking write requests */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR("ncmpi_wait_all")
        end_t = MPI_Wtime();
        timing[3] += end_t - start_t;
    }

    /* tell PnetCDF to release the nonblocking buffer */
    err = ncmpi_buffer_detach(ncid);
    CHECK_ERR("ncmpi_buffer_detach")

    /* obtain the accumulated data amount written by this rank */
    err = ncmpi_inq_put_size(ncid, &w_size);
    CHECK_ERR("ncmpi_inq_put_size")

    /* get all the hints used */
    err = ncmpi_get_file_info(ncid, &info_used);
    CHECK_ERR("ncmpi_get_file_info")

    /* close the output file */
    err = ncmpi_close(ncid);
    CHECK_ERR("ncmpi_close")
    timing[0] = MPI_Wtime() - timing[0];

    /* process timing measurement */
    MPI_Reduce(timing, max_t, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&w_size, &sum_w_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        char value[MPI_MAX_INFO_VAL+1];
        int  flag;

        printf("-----------------------------------------------------------\n");
        printf("---- WRF-IO write benchmark ----\n");
        printf("Output NetCDF file name:        %s\n", out_file);
        printf("Number of MPI processes:        %d\n", nprocs);
        printf("MPI processes arranged to 2D :  %d x %d\n", psizes[0], psizes[1]);
        printf("Grid size longitude x latitude: %lld x %lld\n",longitude,latitude);
        printf("Total number of variables:      %d\n", nvars);
        printf("Number of time records:         %d\n",ntimes);
        printf("Total write amount:             %lld B\n", sum_w_size);
        printf("                                %.2f MiB\n", (float)sum_w_size/1048576);
        printf("                                %.2f GiB\n", (float)sum_w_size/1073741824);
        double bw = (double)sum_w_size / 1048576;
        printf("Max open-to-close time:         %.4f sec\n", max_t[0]);
        printf("Max define metadata time:       %.4f sec\n", max_t[1]);
        printf("Max bput posting  time:         %.4f sec\n", max_t[2]);
        printf("Max wait_all      time:         %.4f sec\n", max_t[3]);
        printf("Write bandwidth:                %.2f MiB/s\n", bw/max_t[0]);
        printf("                                %.2f GiB/s\n", bw/1024.0/max_t[0]);
        printf("-----------------------------------------------------------\n");
        MPI_Info_get(info_used, "striping_factor",  MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_factor:        %s\n", HINT);
        MPI_Info_get(info_used, "striping_unit",    MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_unit:          %s\n", HINT);
        MPI_Info_get(info_used, "cb_buffer_size",   MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_buffer_size:         %s\n", HINT);
        MPI_Info_get(info_used, "cb_nodes",         MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_nodes:               %s\n", HINT);
        MPI_Info_get(info_used, "cb_config_list",   MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_config_list:         %s\n", HINT);
        MPI_Info_get(info_used, "cb_node_list",     MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_node_list:           %s\n", HINT);
        MPI_Info_get(info_used, "nc_pncio",         MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_pncio:              %s\n", HINT);
        MPI_Info_get(info_used, "nc_num_aggrs_per_node",MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_num_aggrs_per_node: %s\n", HINT);
        MPI_Info_get(info_used, "nc_ina_node_list", MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_ina_node_list:      %s\n", HINT);
        MPI_Info_get(info_used, "cray_cb_nodes_multiplier", MPI_MAX_INFO_VAL, value, &flag);
        printf("Hint cray_cb_nodes_multiplier:      %s\n", HINT);
        MPI_Info_get(info_used, "cray_cb_write_lock_mode", MPI_MAX_INFO_VAL, value, &flag);
        printf("Hint cray_cb_write_lock_mode:       %s\n", HINT);
        printf("-----------------------------------------------------------\n");
    }
    MPI_Info_free(&info_used);

err_out:
    if (vars != NULL) {
        for (i=0; i<nvars; i++)
            if (vars[i].nelems > 0)
                free(vars[i].buf);
        free(vars);
    }
    if (err != NC_NOERR) return err;

    if (err != NC_NOERR) return err;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_ENOTENABLED) /* --enable-profiling is not set at configure */
        return NC_NOERR;
    else if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }
    /* report the PnetCDF internal heap memory allocation high water mark */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
        if (verbose && rank == 0)
            printf("Max heap memory allocated by PnetCDF internally is %.2f MiB\n\n",
                   (float)sum_size/1048576);
    }
    fflush(stdout);

    return err;
}

static
int wrf_r_benchmark(char       *in_file,
                    int         psizes[2],
                    int         ntimes,
                    MPI_Info    info)
{
    int i, j, err=NC_NOERR, nprocs, rank, ncid, nvars, dimid;
    double timing[4], max_t[4];
    MPI_Offset buf_size, r_size, sum_r_size, longitude, latitude;
    MPI_Info info_used;
    WRF_VAR *vars=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* open input file */
    err = ncmpi_open(MPI_COMM_WORLD, in_file, NC_NOWRITE, info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line=%d: opening file %s (%s)\n",
               __LINE__, in_file, ncmpi_strerror(err));
        goto err_out;
    }

    /* start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();

    err = ncmpi_inq_dimid(ncid, "south_north", &dimid);
    CHECK_ERR("ncmpi_inq_dimid")
    err = ncmpi_inq_dimlen(ncid, dimid, &longitude);
    CHECK_ERR("ncmpi_inq_dimlen")
    err = ncmpi_inq_dimid(ncid, "west_east", &dimid);
    CHECK_ERR("ncmpi_inq_dimid")
    err = ncmpi_inq_dimlen(ncid, dimid, &latitude);
    CHECK_ERR("ncmpi_inq_dimlen")

    if (debug && rank == 0)
        printf("%s at %d: longitude=%lld latitude=%lld\n",__func__,__LINE__,
               longitude,latitude);

    err = ncmpi_inq_nvars(ncid, &nvars);
    CHECK_ERR("ncmpi_inq_nvars")

    vars = (WRF_VAR*) malloc(sizeof(WRF_VAR) * nvars);

    /* populate variable metadata */
    err = inquire_vars(ncid, vars, psizes, longitude, latitude, &buf_size);
    CHECK_ERR("inquire_vars")

    if (debug) {
        printf("%2d: buf_size %lld\n", rank, buf_size);
        fflush(stdout);
    }

    /* allocate and initialize read buffers */
    MPI_Offset mem_alloc;
    if (debug) mem_alloc = 0;

    for (i=0; i<nvars; i++) {
        if (vars[i].nelems == 0) continue;

        if (vars[i].xtype == NC_FLOAT) {
            float *flt_ptr = (float*) malloc(sizeof(float) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                flt_ptr[j] = rank + j;
            vars[i].buf = (void*) flt_ptr;
            if (debug) mem_alloc += sizeof(float) * vars[i].nelems;
        }
        else if (vars[i].xtype == NC_INT) {
            int *int_ptr = (int*) malloc(sizeof(int) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                int_ptr[j] = rank + j;
            vars[i].buf = (void*) int_ptr;
            if (debug) mem_alloc += sizeof(int) * vars[i].nelems;
        }
        else if (vars[i].xtype == NC_CHAR) {
            char *str_ptr = (char*) malloc(vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                str_ptr[j] = '0' + rank + j;
            vars[i].buf = (void*) str_ptr;
            if (debug) mem_alloc += vars[i].nelems;
        }
    }

    if (debug) {
        if (mem_alloc != buf_size)
            printf("%2d: mem_alloc %lld != buf_size %lld\n", rank, mem_alloc,
                   buf_size);
        fflush(stdout);
        assert(mem_alloc == buf_size);
    }

    timing[1] = MPI_Wtime() - timing[0];
    timing[2] = timing[3] = 0;

    /* ntimes is the number of records */
    for (j=0; j<ntimes; j++) {
        double start_t, end_t;
        start_t = MPI_Wtime();

        if (debug && rank == 0) {
            printf("Reading record %d\n",j);
            fflush(stdout);
        }

        for (i=0; i<nvars; i++) {
            if (vars[i].nelems == 0) continue;

            /* set record ID */
/* TODO: check number of records in input file, must be >= ntimes */
            vars[i].start[0] = j;

            if (vars[i].xtype == NC_FLOAT)
                err = ncmpi_iget_vara_float(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].xtype == NC_INT)
                err = ncmpi_iget_vara_int(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].xtype == NC_CHAR)
                err = ncmpi_iget_vara_text(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            CHECK_ERR(vars[i].name)

        }
        end_t = MPI_Wtime();
        timing[2] += end_t - start_t;
        start_t = end_t;

        if (debug && rank == 0) {
            printf("Flush read requests at end of iteration j=%d\n",j);
            fflush(stdout);
        }

        /* flush all nonblocking read requests */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR("ncmpi_wait_all")
        end_t = MPI_Wtime();
        timing[3] += end_t - start_t;
    }

    /* obtain the accumulated data amount read by this rank */
    err = ncmpi_inq_get_size(ncid, &r_size);
    CHECK_ERR("ncmpi_inq_get_size")

    /* get all the hints used */
    err = ncmpi_get_file_info(ncid, &info_used);
    CHECK_ERR("ncmpi_get_file_info")

    /* close the output file */
    err = ncmpi_close(ncid);
    CHECK_ERR("ncmpi_close")
    timing[0] = MPI_Wtime() - timing[0];

    /* process timing measurement */
    MPI_Reduce(timing, max_t, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r_size, &sum_r_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        char value[MPI_MAX_INFO_VAL+1];
        int  flag;

        printf("-----------------------------------------------------------\n");
        printf("---- WRF-IO read benchmark ----\n");
        printf("Input NetCDF file name:         %s\n", in_file);
        printf("Number of MPI processes:        %d\n", nprocs);
        printf("MPI processes arranged to 2D:   %d x %d\n", psizes[0], psizes[1]);
        printf("Grid size longitude x latitude: %lld x %lld\n",longitude,latitude);
        printf("Total number of variables:      %d\n", nvars);
        printf("Number of time records:         %d\n",ntimes);
        printf("Total read  amount:             %lld B\n", sum_r_size);
        printf("                                %.2f MiB\n", (float)sum_r_size/1048576);
        printf("                                %.2f GiB\n", (float)sum_r_size/1073741824);
        double bw = (double)sum_r_size / 1048576;
        printf("Max open-to-close time:         %.4f sec\n", max_t[0]);
        printf("Max inquire metadata time:      %.4f sec\n", max_t[1]);
        printf("Max iget posting  time:         %.4f sec\n", max_t[2]);
        printf("Max wait_all      time:         %.4f sec\n", max_t[3]);
        printf("Read  bandwidth:                %.2f MiB/s\n", bw/max_t[0]);
        printf("                                %.2f GiB/s\n", bw/1024.0/max_t[0]);
        printf("-----------------------------------------------------------\n");
        MPI_Info_get(info_used, "striping_factor",  MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_factor:        %s\n", HINT);
        MPI_Info_get(info_used, "striping_unit",    MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_unit:          %s\n", HINT);
        MPI_Info_get(info_used, "cb_buffer_size",   MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_buffer_size:         %s\n", HINT);
        MPI_Info_get(info_used, "cb_nodes",         MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_nodes:               %s\n", HINT);
        MPI_Info_get(info_used, "cb_config_list",   MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_config_list:         %s\n", HINT);
        MPI_Info_get(info_used, "cb_node_list",     MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_node_list:           %s\n", HINT);
        MPI_Info_get(info_used, "nc_pncio",         MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_pncio:              %s\n", HINT);
        MPI_Info_get(info_used, "nc_num_aggrs_per_node",MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_num_aggrs_per_node: %s\n", HINT);
        MPI_Info_get(info_used, "nc_ina_node_list", MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_ina_node_list:      %s\n", HINT);
        MPI_Info_get(info_used, "cray_cb_nodes_multiplier", MPI_MAX_INFO_VAL, value, &flag);
        printf("Hint cray_cb_nodes_multiplier:      %s\n", HINT);
        MPI_Info_get(info_used, "cray_cb_write_lock_mode", MPI_MAX_INFO_VAL, value, &flag);
        printf("Hint cray_cb_write_lock_mode:       %s\n", HINT);
        printf("-----------------------------------------------------------\n");
    }
    MPI_Info_free(&info_used);

err_out:
    if (vars != NULL) {
        for (i=0; i<nvars; i++) {
            if (vars[i].nelems > 0) free(vars[i].buf);
            if (vars[i].name != NULL) free(vars[i].name);
            if (vars[i].dimids != NULL) free(vars[i].dimids);
        }
        free(vars);
    }
    if (err != NC_NOERR) return err;

    if (err != NC_NOERR) return err;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_ENOTENABLED) /* --enable-profiling is not set at configure */
        return NC_NOERR;
    else if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }
    /* report the PnetCDF internal heap memory allocation high water mark */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
        if (verbose && rank == 0)
            printf("Max heap memory allocated by PnetCDF internally is %.2f MiB\n\n",
                   (float)sum_size/1048576);
    }
    fflush(stdout);

    return err;
}

/*---- parse_str() >---------------------------------------------------------*/
/* This subroutine parses an input string, in_str, into substring tokens,
 * separated by comma, and returns the number of substrings.
 */
static
int parse_str(char   *in_str,
              char ***str_arr)
{
    char *token, *str_dup;
    int nelems=0;

    str_dup = strdup(in_str);
    if (str_dup == NULL) return nelems;

    /* first find out how many elements there are in in_str */
    token = strtok(str_dup, ",");
    if (token == NULL) {
        free(str_dup);
        return nelems;
    }

    nelems = 1;
    while ((token = strtok(NULL, ",")) != NULL)
        nelems++;

    free(str_dup);

    /* allocate str_arr */
    *str_arr = (char**) malloc(sizeof(char*) * nelems);

    /* populate str_arr[] */
    str_dup = strdup(in_str);
    token = strtok(str_dup, ",");
    (*str_arr)[0] = strdup(token);
    nelems = 1;
    while ((token = strtok(NULL, ",")) != NULL)
        (*str_arr)[nelems++] = strdup(token);

    free(str_dup);
    return nelems;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]\n"
    "       [-h] print this help\n"
    "       [-q] quiet mode\n"
    "       [-d] debug mode\n"
    "       [-r file1,file2,...] benchmark read  performance\n"
    "       [-w file1,file2,...] benchmark write performance\n"
    "       [-y num] longitude of global 2D grid\n"
    "       [-x num] latitude of global 2D grid\n"
    "       [-n num] number of time steps\n"
    "       [-i cdf_file] input text file containing CDL header\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char *out_files, *in_files, *cdl_file, **fname;
    int i, err, nerrs=0, nprocs, rank, ntimes, psizes[2], hid;
    int nfiles, do_read, do_write;
    MPI_Offset longitude, latitude;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose   = 1;
    debug     = 0;
    do_read   = 0;
    do_write  = 0;
    ntimes    = 1;
    cdl_file  = NULL;
    in_files  = NULL;
    out_files = NULL;
    longitude = -1; /* default to use west_east from cdl file */
    latitude  = -1; /* default to use south_north from cdl file */

    while ((i = getopt(argc, argv, "hqdr:w:y:x:n:c:i:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'd': debug = 1;
                      break;
            case 'r': do_read = 1;
                      in_files = strdup(optarg);
                      break;
            case 'w': do_write = 1;
                      out_files = strdup(optarg);
                      break;
            case 'y': longitude = atoll(optarg);
                      break;
            case 'x': latitude = atoll(optarg);
                      break;
            case 'n': ntimes = atoi(optarg);
                      break;
            case 'i': cdl_file = strdup(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    /* check read or write benchmark */
    if (do_read == 0 && do_write == 0) {
        if (rank == 0) {
            fprintf(stderr, "Error: must select read or write benchmark by setting -r and/or -w\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    /* input CDL file is required for write benchmark */
    if (do_write && cdl_file == NULL) {
        if (rank == 0) {
            fprintf(stderr, "Error: write benchmark requires input CDL file\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    /* set up the 2D block-block data partitioning pattern */
    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    if (do_write) {
        /* read and parse the input CDL header file */
        err = cdl_hdr_open(cdl_file, &hid);
        free(cdl_file);
        if (err != NC_NOERR) goto err_out;

        if (debug && rank == 0) {
            printf("longitude=%lld latitude=%lld psizes=%d x %d\n",
                longitude,latitude,psizes[0],psizes[1]);
            fflush(stdout);
        }

        /* Example of out_files: "0.nc,1.nc,2.nc", i.e. 3 output files */
        nfiles = parse_str(out_files, &fname);

        for (i=0; i<nfiles; i++) {
            MPI_Barrier(MPI_COMM_WORLD);

            err = wrf_w_benchmark(fname[i], hid, psizes, longitude, latitude,
                                  ntimes, info);

            if (err != NC_NOERR) goto err_out;
            free(fname[i]);
        }

        if (nfiles > 0) free(fname);

        /* close the CDL file */
        cdl_hdr_close(hid);
        if (out_files != NULL) free(out_files);
    }

    if (do_read) {
        /* Example of in_files: "0.nc,1.nc,2.nc", i.e. 3 input files */
        nfiles = parse_str(in_files, &fname);

        for (i=0; i<nfiles; i++) {
            MPI_Barrier(MPI_COMM_WORLD);

            err = wrf_r_benchmark(fname[i], psizes, ntimes, info);

            if (err != NC_NOERR) goto err_out;
            free(fname[i]);
        }

        if (nfiles > 0) free(fname);
        if (in_files  != NULL) free(in_files);
    }

err_out:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    MPI_Finalize();
    return (nerrs > 0);
}

