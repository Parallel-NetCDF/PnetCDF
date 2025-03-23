/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program evaluates the file write performance of WRF (Wether Research
 * and Forecast Model, https://github.com/wrf-model/WRF) developed at NCAR.
 * It's data partitioning pattern is a 2D block-block checkerboard pattern,
 * along the longitude and latitude. This benchmark program reads a CDL header
 * file in text format (i.e. a netCDF file obtained from "ncmpidump -h") and
 * creates a new file with the same metadata.
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
 *    % mpiexec -n 4 ./wrf_io -l 100 -w 100 -n 2 -i wrf_header.txt ./wrf_io.nc
 *
 *    -----------------------------------------------------------
 *    ---- WRF-IO write benchmark ----
 *    Output NetCDF file name:       ./wrf_io.nc
 *    Number of MPI processes:       4
 *    MPI processes arranged to 2D:  2 x 2
 *    Grid size logitute x latitute: 100 x 100
 *    Total number of variables:     202
 *    Number of time records:        2
 *    Total write amount:            68696794 B
 *                                   65.51 MiB
 *                                   0.06 GiB
 *    Max open-to-close time:        0.0978 sec
 *    Max define metadata time:      0.0085 sec
 *    Max bput posting  time:        0.0148 sec
 *    Max wait_all      time:        0.0718 sec
 *    Write bandwidth:               670.08 MiB/s
 *                                   0.65 GiB/s
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
                   MPI_Offset *logitute,
                   MPI_Offset *latitute,
                   int         psizes[2],
                   MPI_Offset *buf_size)
{
    char *name;
    int i, err=NC_NOERR, nprocs, rank;
    int ndims, *dimids, nvars, my_rank_y, my_rank_x;
    MPI_Offset dimlen, my_start_y, my_start_x, my_count_y, my_count_x;
    nc_type xtype;

    /* determine whether to use logitute, latitute set at command line */
    err = cdl_hdr_inq_ndims(hid, &ndims);
    CHECK_ERR("cdl_hdr_inq_ndims")

    if (*logitute <= 0) {
        for (i=0; i<ndims; i++) {
            err = cdl_hdr_inq_dim(hid, i, &name, &dimlen);
            CHECK_ERR("cdl_hdr_inq_dim")

            /* Note south_north must be defined before south_north_stag */
            if (!strcmp(name, "south_north")) {
                *logitute = dimlen;
                break;
            }
        }
    }
    if (*latitute <= 0) {
        for (i=0; i<ndims; i++) {
            err = cdl_hdr_inq_dim(hid, i, &name, &dimlen);
            CHECK_ERR("cdl_hdr_inq_dim")

            /* Note west_east must be defined before west_east_stag */
            if (!strcmp(name, "west_east")) {
                *latitute = dimlen;
                break;
            }
        }
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    my_rank_y = rank / psizes[1];
    my_rank_x = rank % psizes[1];

    my_count_y = *logitute / psizes[0];
    my_start_y = my_count_y * my_rank_y;
    if (my_rank_y < *logitute % psizes[0]) {
        my_start_y += my_rank_y;
        my_count_y++;
    }
    else {
        my_start_y += *logitute % psizes[0];
    }

    my_count_x = *latitute / psizes[1];
    my_start_x = my_count_x * my_rank_x;
    if (my_rank_x < *latitute % psizes[1]) {
        my_start_x += my_rank_x;
        my_count_x++;
    }
    else {
        my_start_x += *latitute % psizes[1];
    }

    if (debug) {
        printf("%2d: rank (%2d, %2d) start %4lld %4lld count %4lld %4lld\n",
               rank, my_rank_y, my_rank_x, my_start_y, my_start_x, my_count_y, my_count_x);
        fflush(stdout);
    }

    *buf_size = 0;

    /* variables */
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")
    if (debug) printf("var: nvars %d\n", nvars);

    for (i=0; i<nvars; i++, vars++) {
        char *dname1, *dname2, *dname3;
        MPI_Offset dim1, dim2, dim3;

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
int def_dims_vars(int         ncid,
                  int         hid, /* CDL header ID */
                  MPI_Offset *logitute,
                  MPI_Offset *latitute,
                  WRF_VAR    *vars)
{
    char *name;
    void *value;
    int i, j, err=NC_NOERR, ndims, nvars, nattrs;
    MPI_Offset size, nelems;
    nc_type xtype;

    /* define dimensions */
    err = cdl_hdr_inq_ndims(hid, &ndims);
    CHECK_ERR("cdl_hdr_inq_ndims")
    if (debug) printf("dim: ndims %d\n", ndims);

    for (i=0; i<ndims; i++) {
        int dimid;

        err = cdl_hdr_inq_dim(hid, i, &name, &size);
        CHECK_ERR("cdl_hdr_inq_dim")
        if (debug) printf("\t name %s size %lld\n",name, size);

             if (!strcmp(name, "south_north"))      size = *logitute;
        else if (!strcmp(name, "south_north_stag")) size = *logitute + 1;
        else if (!strcmp(name, "west_east"))        size = *latitute;
        else if (!strcmp(name, "west_east_stag"))   size = *latitute + 1;

        err = ncmpi_def_dim(ncid, name, size, &dimid);
        CHECK_ERR("ncmpi_def_dim")
    }

    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")

    /* define variables */
    for (i=0; i<nvars; i++, vars++) {
        err = ncmpi_def_var(ncid, vars->name, vars->xtype, vars->ndims,
                            vars->dimids, &vars->varid);
        CHECK_ERR(name)

        /* define local attributes */
        err = cdl_hdr_inq_nattrs(hid, i, &nattrs);
        CHECK_ERR("cdl_hdr_inq_nattrs")

        if (debug) {
            printf("\t name %s type %d ndims %d nattr %d\n",
                          name, vars->xtype, ndims, nattrs);
            for (j=0; j<ndims; j++)
                printf("\t\tdimid %d\n",vars->dimids[j]);
        }

        for (j=0; j<nattrs; j++) {
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
    err = cdl_hdr_inq_nattrs(hid, NC_GLOBAL, &nattrs);
    CHECK_ERR("cdl_hdr_inq_nattrs")
    if (debug) printf("global attrs: nattrs %d\n", nattrs);

    for (i=0; i<nattrs; i++) {
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
int wrf_io_benchmark(char       *out_file,
                     int         hid,  /* CDL header ID */
                     int         psizes[2],
                     MPI_Offset  logitute,
                     MPI_Offset  latitute,
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
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR("cdl_hdr_inq_nvars")

    vars = (WRF_VAR*) malloc(sizeof(WRF_VAR) * nvars);

    /* populate variable metadata */
    err = construct_vars(hid, vars, &logitute, &latitute, psizes, &buf_size);
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
    CHECK_ERR("ncmpi_create")

    /* define dimension, variables, and attributes */
    err = def_dims_vars(ncid, hid, &logitute, &latitute, vars);
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
            printf("Flush write requests iteration j=%d\n",j);
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
        printf("Output NetCDF file name:       %s\n", out_file);
        printf("Number of MPI processes:       %d\n", nprocs);
        printf("MPI processes arranged to 2D:  %d x %d\n", psizes[0], psizes[1]);
        printf("Grid size logitute x latitute: %lld x %lld\n",logitute,latitute);
        printf("Total number of variables:     %d\n", nvars);
        printf("Number of time records:        %d\n",ntimes);
        printf("Total write amount:            %lld B\n", sum_w_size);
        printf("                               %.2f MiB\n", (float)sum_w_size/1048576);
        printf("                               %.2f GiB\n", (float)sum_w_size/1073741824);
        double bw = (double)sum_w_size / 1048576;
        printf("Max open-to-close time:        %.4f sec\n", max_t[0]);
        printf("Max define metadata time:      %.4f sec\n", max_t[1]);
        printf("Max bput posting  time:        %.4f sec\n", max_t[2]);
        printf("Max wait_all      time:        %.4f sec\n", max_t[3]);
        printf("Write bandwidth:               %.2f MiB/s\n", bw/max_t[0]);
        printf("                               %.2f GiB/s\n", bw/1024.0/max_t[0]);
        printf("-----------------------------------------------------------\n");
        MPI_Info_get(info_used, "striping_factor",  MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_factor:        %s\n", value);
        MPI_Info_get(info_used, "striping_unit",    MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint striping_unit:          %s\n", value);
        MPI_Info_get(info_used, "cb_buffer_size",   MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_buffer_size:         %s\n", value);
        MPI_Info_get(info_used, "aggr_list",        MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint aggr_list:              %s\n", value);
        MPI_Info_get(info_used, "cb_nodes",         MPI_MAX_INFO_VAL, value, &flag);
        printf("MPI-IO hint cb_nodes:               %s\n", value);
        MPI_Info_get(info_used, "nc_num_aggrs_per_node",MPI_MAX_INFO_VAL, value, &flag);
        printf("PnetCDF hint nc_num_aggrs_per_node: %s\n", value);
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
int parse_str(char  *in_str,
              int  **int_arr)
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

    /* allocate int_arr */
    *int_arr = (int*) malloc(sizeof(int) * nelems);

    /* populate int_arr[] */
    str_dup = strdup(in_str);
    token = strtok(str_dup, ",");
    (*int_arr)[0] = atoi(token);
    nelems = 1;
    while ((token = strtok(NULL, ",")) != NULL)
        (*int_arr)[nelems++] = atoi(token);

    free(str_dup);
    return nelems;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS] -i cdf_file output_file\n"
    "       [-h] print this help\n"
    "       [-q] quiet mode\n"
    "       [-d] debug mode\n"
    "       [-l num] logitute of global 2D grid\n"
    "       [-w num] latitute of global 2D grid\n"
    "       [-n num] number of time steps\n"
    "       [-r str] a list of cb_nodes separated by commas\n"
    "       -i cdf_file: input text file containing CDL header \n"
    "       output_file: output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char out_file[1024], *cdl_file, *cb_nodes_str;
    int i, j, err, nerrs=0, nprocs, rank, ntimes, psizes[2], hid;
    int num_cb_nodes, num_intra_nodes, *cb_nodes;
    int nc_num_aggrs_per_node[]={0};
    MPI_Offset logitute, latitute;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose= 1;
    debug  = 0;
    ntimes = 1;
    cdl_file = NULL;
    logitute = -1; /* default to use west_east from cdl file */
    latitute = -1; /* default to use south_north from cdl file */
    cb_nodes     = NULL;
    cb_nodes_str = NULL;

    while ((i = getopt(argc, argv, "hqdl:w:n:r:i:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'd': debug = 1;
                      break;
            case 'l': logitute = atoll(optarg);
                      break;
            case 'w': latitute = atoll(optarg);
                      break;
            case 'r': cb_nodes_str = strdup(optarg);
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
    if (argv[optind] == NULL) strcpy(out_file, "testfile.nc");
    else                      snprintf(out_file, 1024, "%s", argv[optind]);

    /* input CDL file is required */
    if (cdl_file == NULL) {
        if (rank == 0) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* parse CDL header file */
    err = cdl_hdr_open(cdl_file, &hid);
    free(cdl_file);
    if (err != NC_NOERR) goto err_out;

    /* set up the 2D block-block data partitioning pattern */
    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    if (debug && rank == 0) {
        printf("logitute=%lld latitute=%lld psizes=%d x %d\n",
               logitute,latitute,psizes[0],psizes[1]);
        fflush(stdout);
    }

    if (cb_nodes_str != NULL) {
        num_cb_nodes = parse_str(cb_nodes_str, &cb_nodes);
        free(cb_nodes_str);
    }
    else
        num_cb_nodes = 1;

    num_intra_nodes = sizeof(nc_num_aggrs_per_node) / sizeof(int);

    /* set PnetCDF I/O hints */
    MPI_Info_create(&info);

    for (i=0; i<num_cb_nodes; i++) {
        char str[16];

        /* set hint cb_nodes */
        if (cb_nodes != NULL) {
            sprintf(str, "%d", cb_nodes[i]);
            MPI_Info_set(info, "cb_nodes", str);
        }

        for (j=0; j<num_intra_nodes; j++) {
            sprintf(str, "%d", nc_num_aggrs_per_node[j]);
            MPI_Info_set(info, "nc_num_aggrs_per_node", str);

            if (debug && rank == 0) {
                printf("Info cb_nodes set to %d nc_num_aggrs_per_node set to=%d\n",
                       cb_nodes[i],nc_num_aggrs_per_node[j]);
                fflush(stdout);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            err = wrf_io_benchmark(out_file, hid, psizes, logitute, latitute,
                                   ntimes, info);

            if (err != NC_NOERR)
                printf("%d: Error at %s line=%d: i=%d j=%d error=%s\n",
                       rank, argv[0], __LINE__, i, j, ncmpi_strerror(err));
        }
    }
    MPI_Info_free(&info);
    if (cb_nodes != NULL) free(cb_nodes);

err_out:
    cdl_hdr_close(hid);

    MPI_Finalize();
    return (nerrs > 0);
}

