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
 * along the longitude and latitude. This benchmark program writes to 202
 * variables. Among them, 30 are 1D variables. 25 are 2D variables, 122 are 3D,
 * and 25 are 4D. Only 3D and 4D variables are partitioned. 1D and 2D variables
 * are only written by process of rank 0.
 *
 * The grid size can be set by command-line options of -l and -w.
 *
 * To compile:
 *    % mpicc -O2 wrf_io.c -o wrf_io -lpnetcdf
 *
 * An example of run command:
 *    % mpiexec -n 4 ./wrf_io -l 100 -w 100 -n 2 output_file
 *
 *    -----------------------------------------------------------
 *    ---- WRF-IO write benchmark ----
 *    Output NetCDF file name:   testfile.nc
 *    Number of MPI processes:   4
 *    Total number of variables: 202
 *    Grid size:                 100 x 100
 *    Number of time records:    2
 *    Total write amount:        17213594 B
 *                               16.42 MiB
 *                               0.02 GiB
 *    Max open-to-close time:    7.4215 sec
 *    Max bput posting  time:    0.4369 sec
 *    Max wait_all      time:    6.5418 sec
 *    Write bandwidth:           2.21 MiB/s
 *                               0.00 GiB/s
 *    -----------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy(), strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

static int verbose, debug;

#define NVARS 202

#define CHECK_ERR(name) {                               \
    if (err != NC_NOERR) {                              \
        printf("Error at line=%d: name=%s error=%s\n",  \
               __LINE__, name, ncmpi_strerror(err));    \
        goto err_out;                                   \
    }                                                   \
}

#define A(a,b) a##b

#define DEF_DIM(name, val) {                                                  \
    err = ncmpi_def_dim(ncid, name, val, & dim_##val);                        \
    CHECK_ERR(name)                                                           \
}
#define PUT_GATTR_TXT(name, val) {                                            \
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, name, strlen(val), val);        \
    CHECK_ERR(name)                                                           \
}
#define PUT_GATTR_INT(name, val) {                                            \
    int buf = val;                                                            \
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &buf);          \
    CHECK_ERR(name)                                                           \
}
#define PUT_GATTR_FLOAT(name, val) {                                          \
    float buf = val;                                                          \
    err = ncmpi_put_att_float(ncid, NC_GLOBAL, name, NC_FLOAT, 1, &buf);      \
    CHECK_ERR(name)                                                           \
}
#define PUT_ATTR_TXT(name, val) {                                             \
    err = ncmpi_put_att_text(ncid, vars->varid, name, strlen(val), val);      \
    CHECK_ERR(name)                                                           \
}
#define PUT_ATTR_INT1(name, val) {                                            \
    int buf = val;                                                            \
    err = ncmpi_put_att_int(ncid, vars->varid, name, NC_INT, 1, &buf);        \
    CHECK_ERR(name)                                                           \
}
#define DEF_VAR_1D(name, xtype, dim0) {                                       \
    vars++;                                                                   \
    vars->ndims = 1;                                                          \
    vars->var_name = name;                                                    \
    vars->var_xtype = xtype;                                                  \
    vars->nelems = (rank == 0) ? 1 : 0;                                       \
    vars->start[0] = 0;  /* time dimension dim0 is not partitioned */         \
    vars->count[0] = 1;  /* time dimension dim0 is not partitioned */         \
    err = ncmpi_def_var(ncid, name, xtype, 1, &dim_##dim0, &vars->varid);     \
    CHECK_ERR(name)                                                           \
    if (xtype == NC_FLOAT)                                                    \
        *buf_size += sizeof(float) * vars->nelems;                            \
    else if (xtype == NC_INT)                                                 \
        *buf_size += sizeof(int) * vars->nelems;                              \
    else if (xtype == NC_CHAR)                                                \
        *buf_size += vars->nelems;                                            \
}
#define DEF_VAR_2D(name, xtype, dim0, dim1) {                                 \
    int dimids[2] = {dim_##dim0, dim_##dim1};                                 \
    vars++;                                                                   \
    vars->ndims = 2;                                                          \
    vars->var_name = name;                                                    \
    vars->var_xtype = xtype;                                                  \
    vars->nelems = (rank == 0) ? dim1 : 0;                                    \
    vars->start[0] = 0;    /* time dimension */                               \
    vars->count[0] = 1;    /* time dimension */                               \
    vars->start[1] = 0;    /* dimension dim1 is not partitioned */            \
    vars->count[1] = dim1; /* dimension dim1 is not partitioned */            \
    err = ncmpi_def_var(ncid, name, xtype, 2, dimids, &vars->varid);          \
    CHECK_ERR(name)                                                           \
    if (xtype == NC_FLOAT)                                                    \
        *buf_size += sizeof(float) * vars->nelems;                            \
    else if (xtype == NC_INT)                                                 \
        *buf_size += sizeof(int) * vars->nelems;                              \
    else if (xtype == NC_CHAR)                                                \
        *buf_size += vars->nelems;                                            \
}
#define DEF_VAR_3D(name, xtype, dim0, dim1, dim2) {                           \
    int dimids[3] = {dim_##dim0, dim_##dim1, dim_##dim2};                     \
    vars++;                                                                   \
    vars->ndims = 3;                                                          \
    vars->var_name = name;                                                    \
    vars->var_xtype = xtype;                                                  \
    vars->nelems = my_count_y * my_count_x;                                   \
    vars->start[0] = 0; /* time dimension */                                  \
    vars->count[0] = 1; /* time dimension */                                  \
    vars->start[1] = my_start_y;                                              \
    vars->count[1] = my_count_y;                                              \
    vars->start[2] = my_start_x;                                              \
    vars->count[2] = my_count_x;                                              \
    err = ncmpi_def_var(ncid, name, xtype, 3, dimids, &vars->varid);          \
    CHECK_ERR(name)                                                           \
    if (!strcmp(#dim1, "south_north_stag") && my_rank_y == psizes[0]-1)       \
        vars->count[1]++;                                                     \
    if (!strcmp(#dim2, "west_east_stag") && my_rank_x == psizes[1]-1)         \
        vars->count[2]++;                                                     \
    vars->nelems = vars->count[1] *  vars->count[2];                          \
    if (xtype == NC_FLOAT)                                                    \
        *buf_size += sizeof(float) * vars->nelems;                            \
    else if (xtype == NC_INT)                                                 \
        *buf_size += sizeof(int) * vars->nelems;                              \
    else if (xtype == NC_CHAR)                                                \
        *buf_size += vars->nelems;                                            \
}
#define DEF_VAR_4D(name, xtype, dim0, dim1, dim2, dim3) {                     \
    int dimids[4] = {dim_##dim0, dim_##dim1, dim_##dim2, dim_##dim3};         \
    vars++;                                                                   \
    vars->ndims = 4;                                                          \
    vars->var_name = name;                                                    \
    vars->var_xtype = xtype;                                                  \
    vars->start[0] = 0;    /* time dimension */                               \
    vars->count[0] = 1;    /* time dimension */                               \
    vars->start[1] = 0;    /* this dimension is not partitioned */            \
    vars->count[1] = dim1; /* this dimension is not partitioned */            \
    vars->start[2] = my_start_y;                                              \
    vars->count[2] = my_count_y;                                              \
    vars->start[3] = my_start_x;                                              \
    vars->count[3] = my_count_x;                                              \
    err = ncmpi_def_var(ncid, name, xtype, 4, dimids, &vars->varid);          \
    CHECK_ERR(name)                                                           \
    if (!strcmp(#dim2, "south_north_stag") && my_rank_y == psizes[0]-1)       \
        vars->count[2]++;                                                     \
    if (!strcmp(#dim3, "west_east_stag") && my_rank_x == psizes[1]-1)         \
        vars->count[3]++;                                                     \
    vars->nelems = dim1 * vars->count[2] *  vars->count[3];                   \
    if (xtype == NC_FLOAT)                                                    \
        *buf_size += sizeof(float) * vars->nelems;                            \
    else if (xtype == NC_INT)                                                 \
        *buf_size += sizeof(int) * vars->nelems;                              \
    else if (xtype == NC_CHAR)                                                \
        *buf_size += vars->nelems;                                            \
}

typedef struct {
    int varid;
    int ndims;
    int var_xtype;
    char *var_name;
    MPI_Offset nelems;
    MPI_Offset start[4];
    MPI_Offset count[4];
    void *buf;
} WRF_VAR;

static
int def_dims_vars(int         ncid,
                  WRF_VAR    *vars,
                  int         logitute,
                  int         latitute,
                  int         psizes[2],
                  MPI_Offset *buf_size)
{
    int err=NC_NOERR, nprocs, rank;
    int dim_time, dim_DateStrLen, dim_west_east, dim_south_north;
    int dim_bottom_top, dim_bottom_top_stag, dim_soil_layers_stag;
    int dim_west_east_stag, dim_south_north_stag, dim_seed_dim_stag;
    int my_rank_y, my_rank_x;
    MPI_Offset my_start_y, my_start_x, my_count_y, my_count_x;

    MPI_Offset time = NC_UNLIMITED;
    MPI_Offset DateStrLen = 19;
    MPI_Offset west_east = latitute;
    MPI_Offset south_north = logitute;
    MPI_Offset bottom_top = 34;
    MPI_Offset bottom_top_stag = 35;
    MPI_Offset soil_layers_stag = 4;
    MPI_Offset west_east_stag = latitute + 1;
    MPI_Offset south_north_stag = logitute + 1;
    MPI_Offset seed_dim_stag = 8;

    DEF_DIM("Time",             time)
    DEF_DIM("DateStrLen",       DateStrLen);
    DEF_DIM("west_east",        west_east);
    DEF_DIM("south_north",      south_north);
    DEF_DIM("bottom_top",       bottom_top);
    DEF_DIM("bottom_top_stag",  bottom_top_stag);
    DEF_DIM("soil_layers_stag", soil_layers_stag);
    DEF_DIM("west_east_stag",   west_east_stag);
    DEF_DIM("south_north_stag", south_north_stag);
    DEF_DIM("seed_dim_stag",    seed_dim_stag);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    my_rank_y = rank / psizes[1];
    my_rank_x = rank % psizes[1];

    my_count_y = logitute / psizes[0];
    my_start_y = my_count_y * my_rank_y;
    if (my_rank_y < logitute % psizes[0]) {
        my_start_y += my_rank_y;
        my_count_y++;
    }
    else {
        my_start_y += logitute % psizes[0];
    }

    my_count_x = latitute / psizes[1];
    my_start_x = my_count_x * my_rank_x;
    if (my_rank_x < latitute % psizes[1]) {
        my_start_x += my_rank_x;
        my_count_x++;
    }
    else {
        my_start_x += latitute % psizes[1];
    }

    if (debug) {
        printf("%2d: rank (%2d, %2d) start %4lld %4lld count %4lld %4lld\n",
               rank, my_rank_y, my_rank_x, my_start_y, my_start_x, my_count_y, my_count_x);
        fflush(stdout);
    }

    vars--;
    *buf_size = 0;

    DEF_VAR_2D("Times", NC_CHAR, time, DateStrLen)

    DEF_VAR_3D("XLAT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LATITUDE, SOUTH IS NEGATIVE")
    PUT_ATTR_TXT("units","degree_north")
    PUT_ATTR_TXT("stagger","")
    PUT_ATTR_TXT("coordinates","XLONG XLAT")

    DEF_VAR_3D("XLONG", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LONGITUDE, WEST IS NEGATIVE")
    PUT_ATTR_TXT("units","degree_east")
    PUT_ATTR_TXT("stagger","")
    PUT_ATTR_TXT("coordinates","XLONG XLAT")

    DEF_VAR_3D("LU_INDEX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LAND USE CATEGORY")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_2D("ZNU", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "eta values on half (mass) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("ZNW", NC_FLOAT, time, bottom_top_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "eta values on full (w) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ZS", NC_FLOAT, time, soil_layers_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "DEPTHS OF CENTERS OF SOIL LAYERS")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("DZS", NC_FLOAT, time, soil_layers_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "THICKNESSES OF SOIL LAYERS")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_3D("VAR_SSO", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "variance of subgrid-scale orography")
    PUT_ATTR_TXT("units", "m2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("BATHYMETRY_FLAG", NC_INT, time)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "Flag for bathymetry in the global attributes for metgrid data")
    PUT_ATTR_TXT("units", "-")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_4D("U", NC_FLOAT, time, bottom_top, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "x-wind component")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U XTIME")

    DEF_VAR_4D("V", NC_FLOAT, time, bottom_top, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "y-wind component")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V XTIME")

    DEF_VAR_4D("W", NC_FLOAT, time, bottom_top_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "z-wind component")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("PH", NC_FLOAT, time, bottom_top_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "perturbation geopotential")
    PUT_ATTR_TXT("units", "m2 s-2")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("PHB", NC_FLOAT, time, bottom_top_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "base-state geopotential")
    PUT_ATTR_TXT("units", "m2 s-2")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("T", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "perturbation potential temperature theta-t0")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("THM", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "either 1) pert moist pot temp=(1+Rv/Rd Qv)*(theta)-T0, or 2) pert dry pot temp=t")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("HFX_FORCE", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface sensible heat flux")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("LH_FORCE", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface latent heat flux")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("TSK_FORCE", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface skin temperature")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("HFX_FORCE_TEND", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface sensible heat flux tendency")
    PUT_ATTR_TXT("units", "W m-2 s-1")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("LH_FORCE_TEND", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface latent heat flux tendency")
    PUT_ATTR_TXT("units", "W m-2 s-1")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("TSK_FORCE_TEND", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "SCM ideal surface skin temperature tendency")
    PUT_ATTR_TXT("units", "W m-2 s-1")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_3D("MU", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "perturbation dry air mass in column")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("MUB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "base state dry air mass in column")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("NEST_POS", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "-")
    PUT_ATTR_TXT("units", "-")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("P", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "perturbation pressure")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("PB", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "BASE STATE PRESSURE")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_2D("FNM", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "upper weight for vertical stretching")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("FNP", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "lower weight for vertical stretching")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("RDNW", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "inverse d(eta) values between full (w) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("RDN", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "inverse d(eta) values between half (mass) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("DNW", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "d(eta) values between full (w) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("DN", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "d(eta) values between half (mass) levels")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("CFN", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "extrapolation constant")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("CFN1", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "extrapolation constant")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("THIS_IS_AN_IDEAL_RUN", NC_INT, time)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "T/F flag: this is an ARW ideal simulation")
    PUT_ATTR_TXT("units", "-")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_4D("P_HYD", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "hydrostatic pressure")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("Q2", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "QV at 2 M")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("T2", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "TEMP at 2 M")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("TH2", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "POT TEMP at 2 M")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("PSFC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SFC PRESSURE")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("U10", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "U at 10 M")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("V10", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "V at 10 M")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("RDX", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "INVERSE X GRID LENGTH")
    PUT_ATTR_TXT("units", "m-1")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("RDY", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "INVERSE Y GRID LENGTH")
    PUT_ATTR_TXT("units", "m-1")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_3D("AREA2D", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Horizontal grid cell area, using dx, dy, and map factors")
    PUT_ATTR_TXT("units", "m2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("DX2D", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Horizontal grid distance: sqrt(area2d)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("RESM", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "TIME WEIGHT CONSTANT FOR SMALL STEPS")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("ZETATOP", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "ZETA AT MODEL TOP")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("CF1", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "2nd order extrapolation constant")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("CF2", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "2nd order extrapolation constant")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("CF3", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "2nd order extrapolation constant")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("ITIMESTEP", NC_INT, time)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("XTIME", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description","minutes since 2022-06-22 00:00:00")
    PUT_ATTR_TXT("units","minutes since 2022-06-22 00:00:00")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_4D("QVAPOR", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Water vapor mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QCLOUD", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Cloud water mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QRAIN", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Rain water mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QICE", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Ice mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QSNOW", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Snow mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QGRAUP", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Graupel mixing ratio")
    PUT_ATTR_TXT("units", "kg kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QNICE", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Ice Number concentration")
    PUT_ATTR_TXT("units", "  kg-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("QNRAIN", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Rain Number concentration")
    PUT_ATTR_TXT("units", "  kg(-1)")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SHDMAX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ANNUAL MAX VEG FRACTION")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SHDMIN", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ANNUAL MIN VEG FRACTION")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SNOALB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ANNUAL MAX SNOW ALBEDO IN FRACTION")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("TSLB", NC_FLOAT, time, soil_layers_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "SOIL TEMPERATURE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("SMOIS", NC_FLOAT, time, soil_layers_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "SOIL MOISTURE")
    PUT_ATTR_TXT("units", "m3 m-3")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("SH2O", NC_FLOAT, time, soil_layers_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "SOIL LIQUID WATER")
    PUT_ATTR_TXT("units", "m3 m-3")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("SMCREL", NC_FLOAT, time, soil_layers_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "RELATIVE SOIL MOISTURE")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SEAICE", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SEA ICE FLAG")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("XICEM", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SEA ICE FLAG (PREVIOUS STEP)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SFROFF", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SURFACE RUNOFF")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("UDROFF", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "UNDERGROUND RUNOFF")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("IVGTYP", NC_INT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "DOMINANT VEGETATION CATEGORY")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ISLTYP", NC_INT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "DOMINANT SOIL CATEGORY")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("VEGFRA", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "VEGETATION FRACTION")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("GRDFLX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "GROUND HEAT FLUX")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACGRDFLX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED GROUND HEAT FLUX")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSNOM", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED MELTED SNOW")
    PUT_ATTR_TXT("units", "kg m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SNOW", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SNOW WATER EQUIVALENT")
    PUT_ATTR_TXT("units", "kg m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SNOWH", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "PHYSICAL SNOW DEPTH")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("CANWAT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "CANOPY WATER")
    PUT_ATTR_TXT("units", "kg m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SSTSK", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SKIN SEA SURFACE TEMPERATURE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("WATER_DEPTH", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "global water depth")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("COSZEN", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "COS of SOLAR ZENITH ANGLE")
    PUT_ATTR_TXT("units", "dimensionless")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LAI", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LEAF AREA INDEX")
    PUT_ATTR_TXT("units", "m-2/m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("U10E", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Special U at 10 M from MYJSFC")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("V10E", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Special V at 10 M from MYJSFC")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("VAR", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "STANDARD DEVIATION OF SUBGRID-SCALE OROGRAPHY")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("TKE_PBL", NC_FLOAT, time, bottom_top_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "TKE from PBL")
    PUT_ATTR_TXT("units", "m2 s-2")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("EL_PBL", NC_FLOAT, time, bottom_top_stag, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "Length scale from PBL")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "Z")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("MAPFAC_M", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on mass grid")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("MAPFAC_U", NC_FLOAT, time, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on u-grid")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U XTIME")

    DEF_VAR_3D("MAPFAC_V", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on v-grid")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V XTIME")

    DEF_VAR_3D("MAPFAC_MX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on mass grid, x direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("MAPFAC_MY", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on mass grid, y direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("MAPFAC_UX", NC_FLOAT, time, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on u-grid, x direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U XTIME")

    DEF_VAR_3D("MAPFAC_UY", NC_FLOAT, time, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on u-grid, y direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U XTIME")

    DEF_VAR_3D("MAPFAC_VX", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on v-grid, x direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V XTIME")

    DEF_VAR_3D("MF_VX_INV", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Inverse map scale factor on v-grid, x direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V XTIME")

    DEF_VAR_3D("MAPFAC_VY", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Map scale factor on v-grid, y direction")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V XTIME")

    DEF_VAR_3D("F", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Coriolis sine latitude term")
    PUT_ATTR_TXT("units", "s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("E", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Coriolis cosine latitude term")
    PUT_ATTR_TXT("units", "s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SINALPHA", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Local sine of map rotation")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("COSALPHA", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Local cosine of map rotation")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("HGT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "Terrain Height")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("TSK", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SURFACE SKIN TEMPERATURE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("P_TOP", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "PRESSURE TOP OF THE MODEL")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("GOT_VAR_SSO", NC_INT, time)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "whether VAR_SSO was included in WPS output (beginning V3.4)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("T00", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "BASE STATE TEMPERATURE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("P00", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "BASE STATE PRESSURE")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("TLP", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "BASE STATE LAPSE RATE")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("TISO", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "TEMP AT WHICH THE BASE T TURNS CONST")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("TLP_STRAT", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "BASE STATE LAPSE RATE (DT/D(LN(P)) IN STRATOSPHERE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("P_STRAT", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "BASE STATE PRESSURE AT BOTTOM OF STRATOSPHERE")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("MAX_MSFTX", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "Max map factor in domain")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_1D("MAX_MSFTY", NC_FLOAT, time)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "Max map factor in domain")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_3D("RAINC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED TOTAL CUMULUS PRECIPITATION")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("RAINSH", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED SHALLOW CUMULUS PRECIPITATION")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("RAINNC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED TOTAL GRID SCALE PRECIPITATION")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SNOWNC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED TOTAL GRID SCALE SNOW AND ICE")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("GRAUPELNC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED TOTAL GRID SCALE GRAUPEL")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("HAILNC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED TOTAL GRID SCALE HAIL")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_4D("CLDFRA", NC_FLOAT, time, bottom_top, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XYZ")
    PUT_ATTR_TXT("description", "CLOUD FRACTION")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWDOWN", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("GLW", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "DOWNWARD LONG WAVE FLUX AT GROUND SURFACE")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWNORM", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWUPT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWUPTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING CLEAR SKY SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWDNT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWDNTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWUPB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWUPBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWDNB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACSWDNBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWUPT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWUPTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING CLEAR SKY LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWDNT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWDNTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING CLEAR SKY LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWUPB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWUPBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWELLING CLEAR SKY LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWDNB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLWDNBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED DOWNWELLING CLEAR SKY LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWUPT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWUPTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING CLEAR SKY SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWDNT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWDNTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWUPB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWUPBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWDNB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SWDNBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWUPT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWUPTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING CLEAR SKY LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWDNT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWDNTC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING CLEAR SKY LONGWAVE FLUX AT TOP")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWUPB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWUPBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS UPWELLING CLEAR SKY LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWDNB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LWDNBC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "INSTANTANEOUS DOWNWELLING CLEAR SKY LONGWAVE FLUX AT BOTTOM")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("OLR", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "TOA OUTGOING LONG WAVE")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("XLAT_U", NC_FLOAT, time, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LATITUDE, SOUTH IS NEGATIVE")
    PUT_ATTR_TXT("units", "degree_north")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U")

    DEF_VAR_3D("XLONG_U", NC_FLOAT, time, south_north, west_east_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LONGITUDE, WEST IS NEGATIVE")
    PUT_ATTR_TXT("units", "degree_east")
    PUT_ATTR_TXT("stagger", "X")
    PUT_ATTR_TXT("coordinates", "XLONG_U XLAT_U")

    DEF_VAR_3D("XLAT_V", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LATITUDE, SOUTH IS NEGATIVE")
    PUT_ATTR_TXT("units", "degree_north")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V")

    DEF_VAR_3D("XLONG_V", NC_FLOAT, time, south_north_stag, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LONGITUDE, WEST IS NEGATIVE")
    PUT_ATTR_TXT("units", "degree_east")
    PUT_ATTR_TXT("stagger", "Y")
    PUT_ATTR_TXT("coordinates", "XLONG_V XLAT_V")

    DEF_VAR_3D("ALBEDO", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ALBEDO")
    PUT_ATTR_TXT("units", "-")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("CLAT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "COMPUTATIONAL GRID LATITUDE, SOUTH IS NEGATIVE")
    PUT_ATTR_TXT("units", "degree_north")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ALBBCK", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "BACKGROUND ALBEDO")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("EMISS", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SURFACE EMISSIVITY")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("NOAHRES", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "RESIDUAL OF THE NOAH SURFACE ENERGY BUDGET")
    PUT_ATTR_TXT("units", "W m{-2}")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("TMN", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SOIL TEMPERATURE AT LOWER BOUNDARY")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("XLAND", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LAND MASK (1 FOR LAND, 2 FOR WATER)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("UST", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "U* IN SIMILARITY THEORY")
    PUT_ATTR_TXT("units", "m s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("PBLH", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "PBL HEIGHT")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("HFX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "UPWARD HEAT FLUX AT THE SURFACE")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("QFX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "UPWARD MOISTURE FLUX AT THE SURFACE")
    PUT_ATTR_TXT("units", "kg m-2 s-1")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LH", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LATENT HEAT FLUX AT THE SURFACE")
    PUT_ATTR_TXT("units", "W m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACHFX", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWARD HEAT FLUX AT THE SURFACE")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("ACLHF", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "ACCUMULATED UPWARD LATENT HEAT FLUX AT THE SURFACE")
    PUT_ATTR_TXT("units", "J m-2")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SNOWC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SR", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "fraction of frozen precipitation")
    PUT_ATTR_TXT("units", "-")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_1D("SAVE_TOPO_FROM_REAL", NC_INT, time)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "0  ")
    PUT_ATTR_TXT("description", "1=original topo from real/0=topo modified by WRF")
    PUT_ATTR_TXT("units", "flag")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("ISEEDARR_SPPT", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, SPPT")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ISEEDARR_SKEBS", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, SKEBS")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ISEEDARR_RAND_PERTURB", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, RAND_PERT")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ISEEDARRAY_SPP_CONV", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, RAND_PERT2")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ISEEDARRAY_SPP_PBL", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, RAND_PERT3")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("ISEEDARRAY_SPP_LSM", NC_INT, time, seed_dim_stag)
    PUT_ATTR_INT1("FieldType", 106)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "Array to hold seed for restart, RAND_PERT4")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("C1H", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "half levels, c1h = d bf / d eta, using znw")
    PUT_ATTR_TXT("units", "Dimensionless")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("C2H", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "half levels, c2h = (1-c1h)*(p0-pt)")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("C1F", NC_FLOAT, time, bottom_top_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "full levels, c1f = d bf / d eta, using znu")
    PUT_ATTR_TXT("units", "Dimensionless")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("C2F", NC_FLOAT, time, bottom_top_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "full levels, c2f = (1-c1f)*(p0-pt)")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("C3H", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "half levels, c3h = bh")
    PUT_ATTR_TXT("units", "Dimensionless")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("C4H", NC_FLOAT, time, bottom_top)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "half levels, c4h = (eta-bh)*(p0-pt), using znu")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")

    DEF_VAR_2D("C3F", NC_FLOAT, time, bottom_top_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "full levels, c3f = bf")
    PUT_ATTR_TXT("units", "Dimensionless")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_2D("C4F", NC_FLOAT, time, bottom_top_stag)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "Z  ")
    PUT_ATTR_TXT("description", "full levels, c4f = (eta-bf)*(p0-pt), using znw")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "Z")

    DEF_VAR_3D("PCB", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "base state dry air mass in column")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("PC", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "perturbation dry air mass in column")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LANDMASK", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LAND MASK (1 FOR LAND, 0 FOR WATER)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("LAKEMASK", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "LAKE MASK (1 FOR LAKE, 0 FOR NON-LAKE)")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SST", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SEA SURFACE TEMPERATURE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

    DEF_VAR_3D("SST_INPUT", NC_FLOAT, time, south_north, west_east)
    PUT_ATTR_INT1("FieldType", 104)
    PUT_ATTR_TXT("MemoryOrder", "XY ")
    PUT_ATTR_TXT("description", "SEA SURFACE TEMPERATURE FROM WRFLOWINPUT FILE")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("stagger", "")
    PUT_ATTR_TXT("coordinates", "XLONG XLAT XTIME")

// global attributes:
    PUT_GATTR_TXT("TITLE", " OUTPUT FROM WRF V4.4.1 MODEL")
    PUT_GATTR_TXT("START_DATE", "2022-06-22_00:00:00")
    PUT_GATTR_TXT("SIMULATION_START_DATE", "2022-06-22_00:00:00")
    PUT_GATTR_INT("WEST-EAST_GRID_DIMENSION", 1901)
    PUT_GATTR_INT("SOUTH-NORTH_GRID_DIMENSION", 1301)
    PUT_GATTR_INT("BOTTOM-TOP_GRID_DIMENSION", 35)
    PUT_GATTR_FLOAT("DX", 2500.0)
    PUT_GATTR_FLOAT("DY", 2500.0)
    PUT_GATTR_INT("AERCU_OPT", 0)
    PUT_GATTR_FLOAT("AERCU_FCT", 1.0)
    PUT_GATTR_INT("IDEAL_CASE", 0)
    PUT_GATTR_INT("DIFF_6TH_SLOPEOPT", 0)
    PUT_GATTR_INT("AUTO_LEVELS_OPT", 2)
    PUT_GATTR_FLOAT("DIFF_6TH_THRESH", 0.1)
    PUT_GATTR_FLOAT("DZBOT", 50.0)
    PUT_GATTR_FLOAT("DZSTRETCH_S", 1.3)
    PUT_GATTR_FLOAT("DZSTRETCH_U", 1.1)
    PUT_GATTR_INT("SKEBS_ON", 0)
    PUT_GATTR_INT("USE_Q_DIABATIC", 0)
    PUT_GATTR_TXT("GRIDTYPE", "C")
    PUT_GATTR_INT("DIFF_OPT", 1)
    PUT_GATTR_INT("KM_OPT", 4)
    PUT_GATTR_INT("DAMP_OPT", 3)
    PUT_GATTR_FLOAT("DAMPCOEF", 0.2)
    PUT_GATTR_FLOAT("KHDIF", 0.0)
    PUT_GATTR_FLOAT("KVDIF", 0.0)
    PUT_GATTR_INT("MP_PHYSICS", 8)
    PUT_GATTR_INT("RA_LW_PHYSICS", 4)
    PUT_GATTR_INT("RA_SW_PHYSICS", 4)
    PUT_GATTR_INT("SF_SFCLAY_PHYSICS", 2)
    PUT_GATTR_INT("SF_SURFACE_PHYSICS", 2)
    PUT_GATTR_INT("BL_PBL_PHYSICS", 2)
    PUT_GATTR_INT("CU_PHYSICS", 0)
    PUT_GATTR_INT("SF_LAKE_PHYSICS", 0)
    PUT_GATTR_INT("SURFACE_INPUT_SOURCE", 3)
    PUT_GATTR_INT("SST_UPDATE", 0)
    PUT_GATTR_INT("GHG_INPUT", 1)
    PUT_GATTR_INT("GRID_FDDA", 0)
    PUT_GATTR_INT("GFDDA_INTERVAL_M", 0)
    PUT_GATTR_INT("GFDDA_END_H", 0)
    PUT_GATTR_INT("GRID_SFDDA", 0)
    PUT_GATTR_INT("SGFDDA_INTERVAL_M", 0)
    PUT_GATTR_INT("SGFDDA_END_H", 0)
    PUT_GATTR_INT("HYPSOMETRIC_OPT", 2)
    PUT_GATTR_INT("USE_THETA_M", 1)
    PUT_GATTR_INT("GWD_OPT", 0)
    PUT_GATTR_INT("SF_URBAN_PHYSICS", 0)
    PUT_GATTR_INT("SF_SURFACE_MOSAIC", 0)
    PUT_GATTR_INT("SF_OCEAN_PHYSICS", 0)
    PUT_GATTR_INT("SHCU_PHYSICS", 0)
    PUT_GATTR_INT("MFSHCONV", 0)
    PUT_GATTR_INT("FEEDBACK", 1)
    PUT_GATTR_INT("SMOOTH_OPTION", 0)
    PUT_GATTR_FLOAT("SWRAD_SCAT", 1.0)
    PUT_GATTR_INT("W_DAMPING", 1)
    PUT_GATTR_FLOAT("DT", 15.0)
    PUT_GATTR_FLOAT("RADT", 10.0)
    PUT_GATTR_FLOAT("BLDT", 0.0)
    PUT_GATTR_FLOAT("CUDT", 5.0)
    PUT_GATTR_INT("AER_OPT", 0)
    PUT_GATTR_INT("SWINT_OPT", 0)
    PUT_GATTR_INT("AER_TYPE", 1)
    PUT_GATTR_INT("AER_AOD550_OPT", 1)
    PUT_GATTR_INT("AER_ANGEXP_OPT", 1)
    PUT_GATTR_INT("AER_SSA_OPT", 1)
    PUT_GATTR_INT("AER_ASY_OPT", 1)
    PUT_GATTR_FLOAT("AER_AOD550_VAL", 0.12)
    PUT_GATTR_FLOAT("AER_ANGEXP_VAL", 1.3)
    PUT_GATTR_FLOAT("AER_SSA_VAL", 0.85)
    PUT_GATTR_FLOAT("AER_ASY_VAL", 0.9)
    PUT_GATTR_INT("MOIST_ADV_OPT", 1)
    PUT_GATTR_INT("SCALAR_ADV_OPT", 1)
    PUT_GATTR_INT("TKE_ADV_OPT", 1)
    PUT_GATTR_INT("DIFF_6TH_OPT", 0)
    PUT_GATTR_FLOAT("DIFF_6TH_FACTOR", 0.12)
    PUT_GATTR_INT("OBS_NUDGE_OPT", 0)
    PUT_GATTR_FLOAT("BUCKET_MM", -1.0)
    PUT_GATTR_FLOAT("BUCKET_J", -1.0)
    PUT_GATTR_FLOAT("PREC_ACC_DT", 0.0)
    PUT_GATTR_INT("ISFTCFLX", 2)
    PUT_GATTR_INT("ISHALLOW", 0)
    PUT_GATTR_INT("ISFFLX", 1)
    PUT_GATTR_INT("ICLOUD", 1)
    PUT_GATTR_INT("ICLOUD_CU", 0)
    PUT_GATTR_INT("TRACER_PBLMIX", 1)
    PUT_GATTR_INT("SCALAR_PBLMIX", 0)
    PUT_GATTR_INT("YSU_TOPDOWN_PBLMIX", 1)
    PUT_GATTR_INT("GRAV_SETTLING", 0)
    PUT_GATTR_INT("CLDOVRLP", 2)
    PUT_GATTR_INT("IDCOR", 0)
    PUT_GATTR_INT("DFI_OPT", 0)
    PUT_GATTR_INT("NTASKS_X", 32)
    PUT_GATTR_INT("NTASKS_Y", 64)
    PUT_GATTR_INT("NTASKS_TOTAL", 2048)
    PUT_GATTR_TXT("SIMULATION_INITIALIZATION_TYPE", "REAL-DATA CASE")
    PUT_GATTR_INT("WEST-EAST_PATCH_START_UNSTAG", 1)
    PUT_GATTR_INT("WEST-EAST_PATCH_END_UNSTAG", 1900)
    PUT_GATTR_INT("WEST-EAST_PATCH_START_STAG", 1)
    PUT_GATTR_INT("WEST-EAST_PATCH_END_STAG", 1901)
    PUT_GATTR_INT("SOUTH-NORTH_PATCH_START_UNSTAG", 1)
    PUT_GATTR_INT("SOUTH-NORTH_PATCH_END_UNSTAG", 1300)
    PUT_GATTR_INT("SOUTH-NORTH_PATCH_START_STAG", 1)
    PUT_GATTR_INT("SOUTH-NORTH_PATCH_END_STAG", 1301)
    PUT_GATTR_INT("BOTTOM-TOP_PATCH_START_UNSTAG", 1)
    PUT_GATTR_INT("BOTTOM-TOP_PATCH_END_UNSTAG", 34)
    PUT_GATTR_INT("BOTTOM-TOP_PATCH_START_STAG", 1)
    PUT_GATTR_INT("BOTTOM-TOP_PATCH_END_STAG", 35)
    PUT_GATTR_INT("GRID_ID", 1)
    PUT_GATTR_INT("PARENT_ID", 0)
    PUT_GATTR_INT("I_PARENT_START", 1)
    PUT_GATTR_INT("J_PARENT_START", 1)
    PUT_GATTR_INT("PARENT_GRID_RATIO", 1)
    PUT_GATTR_FLOAT("CEN_LAT", 40.00001)
    PUT_GATTR_FLOAT("CEN_LON", -98.0)
    PUT_GATTR_FLOAT("TRUELAT1", 30.0)
    PUT_GATTR_FLOAT("TRUELAT2", 60.0)
    PUT_GATTR_FLOAT("MOAD_CEN_LAT", 40.00001)
    PUT_GATTR_FLOAT("STAND_LON", -98.0)
    PUT_GATTR_FLOAT("POLE_LAT", 90.0)
    PUT_GATTR_FLOAT("POLE_LON", 0.0)
    PUT_GATTR_FLOAT("GMT", 0.0)
    PUT_GATTR_INT("JULYR", 2022)
    PUT_GATTR_INT("JULDAY", 173)
    PUT_GATTR_INT("MAP_PROJ", 1)
    PUT_GATTR_TXT("MAP_PROJ_CHAR", "Lambert Conformal")
    PUT_GATTR_TXT("MMINLU", "MODIFIED_IGBP_MODIS_NOAH")
    PUT_GATTR_INT("NUM_LAND_CAT", 21)
    PUT_GATTR_INT("ISWATER", 17)
    PUT_GATTR_INT("ISLAKE", 21)
    PUT_GATTR_INT("ISICE", 15)
    PUT_GATTR_INT("ISURBAN", 13)
    PUT_GATTR_INT("ISOILWATER", 14)
    PUT_GATTR_INT("HYBRID_OPT", 2)
    PUT_GATTR_FLOAT("ETAC", 0.2)

err_out:
    return err;
}

static
int wrf_io_benchmark(char     *filename,
                     int       psizes[2],
                     int       logitute,
                     int       latitute,
                     int       ntimes,
                     MPI_Info  info)
{
    int i, j, err=NC_NOERR, nprocs, rank;
    int cmode, ncid;
    double timing[3], max_t[3];
    MPI_Offset buf_size, w_size, sum_w_size;
    MPI_Info info_used;
    WRF_VAR *vars=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();

    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR("ncmpi_create")

    vars = (WRF_VAR*) malloc(sizeof(WRF_VAR) * NVARS);

    err = def_dims_vars(ncid, vars, logitute, latitute, psizes, &buf_size);
    CHECK_ERR("def_dims_vars")

    err = ncmpi_enddef(ncid);
    CHECK_ERR("ncmpi_enddef")

    if (debug) {
        printf("%2d: buf_size %lld\n", rank, buf_size);
        fflush(stdout);
    }

    /* allocate and initialize write buffer */
    MPI_Offset mem_alloc;
    if (debug) mem_alloc = 0;

    for (i=0; i<NVARS; i++) {
        if (vars[i].nelems == 0) continue;

        if (vars[i].var_xtype == NC_FLOAT) {
            float *flt_ptr = (float*) malloc(sizeof(float) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                flt_ptr[j] = rank + j;
            vars[i].buf = (void*) flt_ptr;
            if (debug) mem_alloc += sizeof(float) * vars[i].nelems;
        }
        else if (vars[i].var_xtype == NC_INT) {
            int *int_ptr = (int*) malloc(sizeof(int) * vars[i].nelems);
            for (j=0; j<vars[i].nelems; j++)
                int_ptr[j] = rank + j;
            vars[i].buf = (void*) int_ptr;
            if (debug) mem_alloc += sizeof(int) * vars[i].nelems;
        }
        else if (vars[i].var_xtype == NC_CHAR) {
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

    err = ncmpi_buffer_attach(ncid, buf_size);
    CHECK_ERR("ncmpi_buffer_attach")

    timing[1] = timing[2] = 0;

    for (j=0; j<ntimes; j++) {
        double start_t, end_t;
        start_t = MPI_Wtime();

        if (debug && rank == 0) {
            printf("Writing record %d\n",j);
            fflush(stdout);
        }

        for (i=0; i<NVARS; i++) {
            if (vars[i].nelems == 0) continue;

            /* set record ID */
            vars[i].start[0] = j;

            if (vars[i].var_xtype == NC_FLOAT)
                err = ncmpi_bput_vara_float(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].var_xtype == NC_INT)
                err = ncmpi_bput_vara_int(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            else if (vars[i].var_xtype == NC_CHAR)
                err = ncmpi_bput_vara_text(ncid, vars[i].varid, vars[i].start,
                                            vars[i].count, vars[i].buf, NULL);
            CHECK_ERR(vars[i].var_name)

        }
        end_t = MPI_Wtime();
        timing[1] += end_t - start_t;

        if (debug && rank == 0) {
            printf("Flush write requests iteration j=%d\n",j);
            fflush(stdout);
        }

        start_t = end_t;
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR("ncmpi_wait_all")
        end_t = MPI_Wtime();
        timing[2] += end_t - start_t;
    }

    err = ncmpi_buffer_detach(ncid);
    CHECK_ERR("ncmpi_buffer_detach")

    err = ncmpi_inq_put_size(ncid, &w_size);
    CHECK_ERR("ncmpi_inq_put_size")

    /* get all the hints used */
    err = ncmpi_get_file_info(ncid, &info_used);
    CHECK_ERR("ncmpi_get_file_info")

    err = ncmpi_close(ncid); CHECK_ERR("ncmpi_close")
    timing[0] = MPI_Wtime() - timing[0];

    MPI_Reduce(timing, max_t, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&w_size, &sum_w_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        char value[MPI_MAX_INFO_VAL+1];
        int  flag;

        printf("-----------------------------------------------------------\n");
        printf("---- WRF-IO write benchmark ----\n");
        printf("Output NetCDF file name:      %s\n", filename);
        printf("Number of MPI processes:       %d\n", nprocs);
        printf("MPI processes arranged to 2D:  %d x %d\n", psizes[0], psizes[1]);
        printf("Grid size logitute x latitute: %d x %d\n",logitute, latitute);
        printf("Total number of variables:     %d\n", NVARS);
        printf("Number of time records:        %d\n",ntimes);
        printf("Total write amount:            %lld B\n", sum_w_size);
        printf("                               %.2f MiB\n", (float)sum_w_size/1048576);
        printf("                               %.2f GiB\n", (float)sum_w_size/1073741824);
        double bw = (double)sum_w_size / 1048576;
        printf("Max open-to-close time:        %.4f sec\n", max_t[0]);
        printf("Max bput posting  time:        %.4f sec\n", max_t[1]);
        printf("Max wait_all      time:        %.4f sec\n", max_t[2]);
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
        for (i=0; i<NVARS; i++)
            if (vars[i].nelems > 0)
                free(vars[i].buf);
        free(vars);
    }
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
    "Usage: %s [-h|-q|-d|-l logitute |-w latitute|-n num | -r str] file_name\n"
    "       [-h] print this help\n"
    "       [-q] quiet mode\n"
    "       [-d] debug mode\n"
    "       [-l num] logitute of global 2D grid\n"
    "       [-w num] latitute of global 2D grid\n"
    "       [-n num] number of time steps\n"
    "       [-r str] a list of cb_nodes separated by commas\n"
    "       filename: output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[1024], *cb_nodes_str;
    int i, j, err, nerrs=0, nprocs, rank, ntimes, psizes[2];
    int logitute, latitute, num_cb_nodes, num_intra_nodes, *cb_nodes;
    int nc_num_aggrs_per_node[]={0};
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose= 1;
    debug  = 0;
    ntimes = 1;
    logitute = 100;
    latitute = 100;
    cb_nodes     = NULL;
    cb_nodes_str = NULL;

    while ((i = getopt(argc, argv, "hqdl:w:n:r:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'd': debug = 1;
                      break;
            case 'l': logitute = atoi(optarg);
                      break;
            case 'w': latitute = atoi(optarg);
                      break;
            case 'r': cb_nodes_str = strdup(optarg);
                      break;
            case 'n': ntimes = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 1024, "%s", argv[optind]);

    /* set up the 2D block-block data partitioning pattern */
    psizes[0] = psizes[1] = 0;
    MPI_Dims_create(nprocs, 2, psizes);

    if (debug && rank == 0) {
        printf("logitute=%d latitute=%d psizes=%d x %d\n",logitute,latitute,
               psizes[0],psizes[1]);
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

            err = wrf_io_benchmark(filename, psizes, logitute, latitute, ntimes, info);

            if (err != NC_NOERR)
                printf("%d: Error at %s line=%d: i=%d j=%d error=%s\n",
                       rank, argv[0], __LINE__, i, j, ncmpi_strerror(err));
        }
    }
    MPI_Info_free(&info);
    if (cb_nodes != NULL) free(cb_nodes);

    MPI_Finalize();
    return (nerrs > 0);
}

