/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program tests writing variables of NC_CHAR types and re-enter define
 * mode to grow file header size.
 *
 * Remove the line of "#defin PnetCDF" to test with NetCDF4.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define LEN_NAME 101

#define USE_PNETCDF

#ifdef USE_PNETCDF
#include <pnetcdf.h>
#define CREATE(c,f,m,i,r) ncmpi_create(c,f,m,i,r)
#define DEF_DIM           ncmpi_def_dim
#define DEFVAR            ncmpi_def_var
#define ENDDEF            ncmpi_enddef
#define REDEF             ncmpi_redef
#define CLOSE             ncmpi_close
#define STRERROR          ncmpi_strerror
#define PUT_TEXT          ncmpi_put_var_text_all
#define PUT_ATTR_TEXT     ncmpi_put_att_text
#else
#include <netcdf.h>
#include <netcdf_par.h>
#define CREATE(c,f,m,i,r) nc_create_par(f,m,c,i,r)
#define DEF_DIM           nc_def_dim
#define DEFVAR            nc_def_var
#define ENDDEF            nc_enddef
#define REDEF             nc_redef
#define CLOSE             nc_close
#define STRERROR          nc_strerror
#define PUT_TEXT          nc_put_var_text
#define PUT_ATTR_TEXT     nc_put_att_text
#endif

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,STRERROR(err)); \
    } \
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *filename="testfile.nc", *buf, *ptr;
    int err, nerrs=0, rank, nprocs, cmode;
    int ncid, dimids[2], varid;
    int num_el_blk=3, num_side_sets=25;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* create a new file */
    cmode = NC_CLOBBER;
    err  = CREATE(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL,
                  &ncid); CHECK_ERR

    /* define dimensions */
    err = DEF_DIM(ncid, "num_el_blk", num_el_blk, &dimids[0]); CHECK_ERR
    err = DEF_DIM(ncid, "len_name", LEN_NAME, &dimids[1]); CHECK_ERR

    /* define variable eb_names */
    err = DEFVAR(ncid, "eb_names", NC_CHAR, 2, dimids, &varid);
    CHECK_ERR

    err = ENDDEF(ncid); CHECK_ERR

    /* write variable eb_names */
    buf = (char*) calloc(num_el_blk * LEN_NAME, 1);
    ptr = buf;
    strcpy(ptr, "block_1");        ptr += LEN_NAME;
    strcpy(ptr, "block_1_air");    ptr += LEN_NAME;
    strcpy(ptr, "block_1_fluid");
    err = PUT_TEXT(ncid, varid, buf); CHECK_ERR
    free(buf);

    /* add more metedata */
    err = REDEF(ncid); CHECK_ERR

    /* add a global attribute to grow header */
    char *attr = (char*) calloc(512, 1);
    err = PUT_ATTR_TEXT(ncid, NC_GLOBAL, "attr", 512, attr);
    CHECK_ERR
    free(attr);

    /* define variable ss_names */
    err = DEF_DIM(ncid, "num_side_sets", num_side_sets, &dimids[0]);
    CHECK_ERR
    err = DEFVAR(ncid, "ss_names", NC_CHAR, 2, dimids, &varid);
    CHECK_ERR

    err = ENDDEF(ncid); CHECK_ERR

    /* write variable ss_names */
    buf = (char*) calloc(num_side_sets * LEN_NAME, 1);
    ptr = buf;
    strcpy(ptr, "surface_1");        ptr += LEN_NAME;
    strcpy(ptr, "surface_1_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_1_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_2");        ptr += LEN_NAME;
    strcpy(ptr, "surface_2_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_2_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_3");        ptr += LEN_NAME;
    strcpy(ptr, "surface_3_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_3_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_4");        ptr += LEN_NAME;
    strcpy(ptr, "surface_44");       ptr += LEN_NAME;
    strcpy(ptr, "surface_44_air");   ptr += LEN_NAME;
    strcpy(ptr, "surface_44_fluid"); ptr += LEN_NAME;
    strcpy(ptr, "surface_4_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_4_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_5");        ptr += LEN_NAME;
    strcpy(ptr, "surface_5_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_5_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_6");        ptr += LEN_NAME;
    strcpy(ptr, "surface_66");       ptr += LEN_NAME;
    strcpy(ptr, "surface_66_air");   ptr += LEN_NAME;
    strcpy(ptr, "surface_66_fluid"); ptr += LEN_NAME;
    strcpy(ptr, "surface_6_air");    ptr += LEN_NAME;
    strcpy(ptr, "surface_6_fluid");  ptr += LEN_NAME;
    strcpy(ptr, "surface_air_fluid");
    err = PUT_TEXT(ncid, varid, buf); CHECK_ERR

    /* add aother global attribute to grow header */
    err = REDEF(ncid); CHECK_ERR

    char *attr2 = (char*) calloc(512, 1);
    err = PUT_ATTR_TEXT(ncid, NC_GLOBAL, "attr2", 512, attr2);
    CHECK_ERR
    free(attr2);

    err = ENDDEF(ncid); CHECK_ERR

    err = CLOSE(ncid); CHECK_ERR
    free(buf);

    MPI_Finalize();
    return (nerrs > 0);
}

