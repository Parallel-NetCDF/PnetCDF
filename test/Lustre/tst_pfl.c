/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* This program test whether PnetCDF can properly open/create a file in a
 * folder configured with Lustre Progressive File Layout.

First create a new directory and file with PFL (with 5 components):

    mkdir -p /eagle/radix-io/wkliao/pfl
    lfs setstripe -E 4m   -c  4 \
                  -E 1g   -c  8 \
                  -E 10g  -c 16 \
                  -E 100g -c 32 \
                  -E eof  -c -1 \
                  /eagle/radix-io/wkliao/pfl

Copy a netCDF file over. The file should inherit the folder's PFL settings.

    cp /eagle/radix-io/wkliao/map_i_case_fill_1344p.nc \
       /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc

To compile:
    mpicc -Wall -O2 -o tst_pfl ./tst_pfl.c -I/home/wkliao/PnetCDF/Github/include -L/home/wkliao/PnetCDF/Github/lib -lpnetcdf -llustreapi

To run:
    ./tst_pfl /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32

==== test OPEN nc_driver = mpiio ==========================
---- print_info ---------------------------
MPI File Info: nkeys = 39
MPI File Info: [15] key =           striping_factor, value = 32
MPI File Info: [22] key =                  cb_nodes, value = 1
MPI File Info: [32] key =                 nc_driver, value = mpiio
MPI File Info: [34] key =             file_striping, value = auto
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc
lmm_stripe_size   = 1048576
lmm_stripe_count  = 8
lmm_stripe_offset = 0

==== test OPEN nc_driver = gio ============================
---- print_info ---------------------------
MPI File Info: nkeys = 30
MPI File Info: [ 0] key =                 nc_driver, value = gio
MPI File Info: [ 9] key =                  cb_nodes, value = 1
MPI File Info: [12] key =           striping_factor, value = 32
MPI File Info: [15] key =        overstriping_ratio, value = 1
MPI File Info: [26] key =             file_striping, value = auto
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc
lmm_stripe_size   = 1048576
lmm_stripe_count  = 8
lmm_stripe_offset = 0



==== test CREATE nc_driver = mpiio ==========================
---- print_info ---------------------------
MPI File Info: nkeys = 39
MPI File Info: [15] key =           striping_factor, value = 1
MPI File Info: [22] key =                  cb_nodes, value = 1
MPI File Info: [32] key =                 nc_driver, value = mpiio
MPI File Info: [34] key =             file_striping, value = auto
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc.mpiio.noinherit
lmm_stripe_size   = 1048576
lmm_stripe_count  = 1
lmm_stripe_offset = 0

==== test CREATE nc_driver = mpiio ==========================
---- print_info ---------------------------
MPI File Info: nkeys = 39
MPI File Info: [15] key =           striping_factor, value = 32
MPI File Info: [22] key =                  cb_nodes, value = 1
MPI File Info: [32] key =                 nc_driver, value = mpiio
MPI File Info: [34] key =             file_striping, value = inherit
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc.mpiio.inherit
lmm_stripe_size   = 1048576
lmm_stripe_count  = 4
lmm_stripe_offset = 0

==== test CREATE nc_driver = gio ============================
---- print_info ---------------------------
MPI File Info: nkeys = 30
MPI File Info: [ 0] key =                 nc_driver, value = gio
MPI File Info: [ 1] key =             file_striping, value = auto
MPI File Info: [10] key =                  cb_nodes, value = 1
MPI File Info: [13] key =           striping_factor, value = 1
MPI File Info: [16] key =        overstriping_ratio, value = 1
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc.gio.noinherit
lmm_stripe_size   = 1048576
lmm_stripe_count  = 1
lmm_stripe_offset = 0

==== test CREATE nc_driver = gio ============================
---- print_info ---------------------------
MPI File Info: nkeys = 30
MPI File Info: [ 0] key =                 nc_driver, value = gio
MPI File Info: [ 1] key =             file_striping, value = inherit
MPI File Info: [10] key =                  cb_nodes, value = 1
MPI File Info: [13] key =           striping_factor, value = 32
MPI File Info: [16] key =        overstriping_ratio, value = 1
---- get_stripe ---------------------------
filename          = /eagle/radix-io/wkliao/pfl/pfl_4_8_16_32.nc.gio.inherit
lmm_stripe_size   = 1048576
lmm_stripe_count  = 4
lmm_stripe_offset = 0

 */

#include <stdio.h>
#include <string.h>
#include <unistd.h> /* unlink() */

#include <mpi.h>
#include <pnetcdf.h>

/*----< print_info() >------------------------------------------------------*/
static
void print_info(MPI_Info *info_used)
{
    int  i, nkeys;

    printf("---- %s ---------------------------\n",__func__);
    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);

        if (strcmp(key, "striping_factor") &&
            strcmp(key, "cb_nodes") &&
            strcmp(key, "file_striping") &&
            strcmp(key, "overstriping_ratio") &&
            strcmp(key, "nc_driver")) continue;

        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}

#define ERR {if(err!=NC_NOERR){fprintf(stderr,"Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <lustre/lustreapi.h>
#define MAX_LOV_UUID_COUNT      1000

void get_stripe(char *filename)
{
    int lumlen, err, fd;
    struct lov_user_md *lum = NULL;

    lumlen = sizeof(struct lov_user_md) + MAX_LOV_UUID_COUNT * sizeof(struct lov_user_ost_data);
    lum = (struct lov_user_md *) calloc(1, lumlen);

    lum->lmm_magic = LOV_USER_MAGIC;

    fd = open(filename, O_RDONLY, 0666);

    err = ioctl(fd, LL_IOC_LOV_GETSTRIPE, (void *) lum);
    if (!err) {
        printf("---- %s ---------------------------\n",__func__);
        printf("filename          = %s\n", filename);
        printf("lmm_stripe_size   = %u\n", lum->lmm_stripe_size);
        printf("lmm_stripe_count  = %u\n", lum->lmm_stripe_count);
        printf("lmm_stripe_offset = %u\n", lum->lmm_stripe_offset);
    }

    close(fd);
}

int test_open(char *filename, MPI_Info info)
{
    int ncid, err, nerrs=0;
    MPI_Info info_used;

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    ERR

    err = ncmpi_inq_file_info(ncid, &info_used);
    ERR
    print_info(&info_used);
    MPI_Info_free(&info_used);

    err = ncmpi_close(ncid);

    get_stripe(filename);

    return nerrs;
}

int test_create(char *filename, MPI_Info info)
{
    int ncid, err, nerrs=0;
    MPI_Info info_used;

    unlink(filename);

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid);
    ERR

    err = ncmpi_inq_file_info(ncid, &info_used);
    ERR
    print_info(&info_used);
    MPI_Info_free(&info_used);

    err = ncmpi_close(ncid);

    get_stripe(filename);

    return nerrs;
}

int main(int argc, char **argv)
{
    char filename[1024];
    int nerrs=0, rank;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 2) {
        if (rank == 0) printf("Usage: %s filename\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    MPI_Info_create(&info);

    printf("\n==== test OPEN nc_driver = mpiio ==========================\n");
    sprintf(filename, "%s.nc", argv[1]);
    MPI_Info_set(info, "nc_driver", "mpiio");
    nerrs += test_open(filename, info);

    printf("\n==== test OPEN nc_driver = gio ============================\n");
    sprintf(filename, "%s.nc", argv[1]);
    MPI_Info_set(info, "nc_driver", "gio");
    nerrs += test_open(filename, info);

    printf("\n\n");
    printf("\n==== test CREATE nc_driver = mpiio ==========================\n");
    sprintf(filename, "%s.nc.mpiio.noinherit", argv[1]);
    MPI_Info_set(info, "nc_driver", "mpiio");
    nerrs += test_create(filename, info);

    printf("\n==== test CREATE nc_driver = mpiio ==========================\n");
    sprintf(filename, "%s.nc.mpiio.inherit", argv[1]);
    MPI_Info_set(info, "nc_driver", "mpiio");
    MPI_Info_set(info, "file_striping", "inherit");
    nerrs += test_create(filename, info);

    printf("\n==== test CREATE nc_driver = gio ============================\n");
    sprintf(filename, "%s.nc.gio.noinherit", argv[1]);
    MPI_Info_set(info, "nc_driver", "gio");
    MPI_Info_set(info, "file_striping", "auto");
    nerrs += test_create(filename, info);

    printf("\n==== test CREATE nc_driver = gio ============================\n");
    sprintf(filename, "%s.nc.gio.inherit", argv[1]);
    MPI_Info_set(info, "nc_driver", "gio");
    MPI_Info_set(info, "file_striping", "inherit");
    nerrs += test_create(filename, info);

    MPI_Info_free(&info);

    MPI_Finalize();

    return 0;
}

