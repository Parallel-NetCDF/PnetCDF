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
#define DEFVAR_FILL       ncmpi_def_var_fill
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
#define DEFVAR_FILL       nc_def_var_fill
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

#define ADD_DIM(name, len) \
    int name = len, name##_d; \
    err = DEF_DIM(ncid, #name, len, &name##_d); \
    CHECK_ERR

#define ADD_VAR_1D(name, xtype, x) \
    int name##_v; \
    dims[0] = x##_d; \
    err = DEFVAR(ncid, #name, xtype, 1, dims, &name##_v); \
    CHECK_ERR \
    err = DEFVAR_FILL(ncid, name##_v, 0, NULL); \
    CHECK_ERR

#define ADD_VAR_2D(name, xtype, y, x) \
    int name##_v; \
    dims[0] = y##_d; dims[1] = x##_d; \
    err = DEFVAR(ncid, #name, xtype, 2, dims, &name##_v); \
    CHECK_ERR \
    err = DEFVAR_FILL(ncid, name##_v, 0, NULL); \
    CHECK_ERR

#define ADD_VAR_3D(name, xtype, z, y, x) \
    int name##_v; \
    dims[0] = z##_d; dims[1] = y##_d; dims[2] = x##_d; \
    err = DEFVAR(ncid, #name, xtype, 3, dims, &name##_v); \
    CHECK_ERR \
    err = DEFVAR_FILL(ncid, name##_v, 0, NULL); \
    CHECK_ERR

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *filename="testfile.nc", *buf, *ptr;
    int err, nerrs=0, rank, nprocs, cmode;
    int ncid, varid[256], nvars=0, dims[4];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* create a new file */
    cmode = NC_CLOBBER;
    err = CREATE(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* define dimensions */
    ADD_DIM(num_qa_rec, 1);
    ADD_DIM(four, 4);
    ADD_DIM(len_string, 33);
    ADD_DIM(num_info, 433);
    ADD_DIM(len_line, 81);
    ADD_DIM(len_name, 101);
    ADD_DIM(num_dim, 3);
    ADD_DIM(time_step, NC_UNLIMITED);
    ADD_DIM(num_nodes, 148977);
    ADD_DIM(num_elem, 796986);
    ADD_DIM(num_el_blk, 3);
    ADD_DIM(num_side_sets, 25);
    ADD_DIM(num_el_in_blk2, 791921);
    ADD_DIM(num_nod_per_el2, 4);
    ADD_DIM(num_el_in_blk3, 5065);
    ADD_DIM(num_nod_per_el3, 4);
    ADD_DIM(num_side_ss2, 24564);
    ADD_DIM(num_side_ss5, 11888);
    ADD_DIM(num_side_ss9, 72);
    ADD_DIM(num_side_ss12, 156);
    ADD_DIM(num_side_ss14, 779);
    ADD_DIM(num_side_ss15, 515);
    ADD_DIM(num_side_ss17, 24540);
    ADD_DIM(num_side_ss21, 36);
    ADD_DIM(num_side_ss25, 1894);
    ADD_DIM(num_elem_maps, 1);

    /* define variables */
    ADD_VAR_3D(qa_records, NC_CHAR, num_qa_rec, four, len_string);
    ADD_VAR_2D(info_records,NC_CHAR, num_info, len_line) ;
    ADD_VAR_1D(time_whole, NC_DOUBLE, time_step);
    ADD_VAR_1D(node_num_map, NC_INT, num_nodes) ;
    ADD_VAR_1D(elem_num_map, NC_INT, num_elem) ;
    ADD_VAR_1D(eb_status, NC_INT, num_el_blk) ;
    ADD_VAR_1D(eb_prop1, NC_INT, num_el_blk) ;
    err = PUT_ATTR_TEXT(ncid, eb_prop1_v, "name", 2, "ID"); CHECK_ERR
    ADD_VAR_2D(eb_names, NC_CHAR, num_el_blk, len_name);
    ADD_VAR_1D(ss_status, NC_INT, num_side_sets) ;
    ADD_VAR_1D(ss_prop1, NC_INT, num_side_sets) ;
    err = PUT_ATTR_TEXT(ncid, ss_prop1_v, "name", 2, "ID"); CHECK_ERR
    ADD_VAR_2D(ss_names, NC_CHAR, num_side_sets, len_name) ;
    ADD_VAR_1D(coordx, NC_DOUBLE, num_nodes);
    ADD_VAR_1D(coordy, NC_DOUBLE, num_nodes);
    ADD_VAR_1D(coordz, NC_DOUBLE, num_nodes);
    ADD_VAR_2D(coor_names, NC_CHAR, num_dim, len_name) ;
    ADD_VAR_2D(connect2, NC_INT, num_el_in_blk2, num_nod_per_el2) ;
    err = PUT_ATTR_TEXT(ncid, connect2_v, "elem_type", 6, "tetra4"); CHECK_ERR
    ADD_VAR_2D(connect3, NC_INT, num_el_in_blk3, num_nod_per_el3) ;
    err = PUT_ATTR_TEXT(ncid, connect3_v, "elem_type", 6, "tetra4"); CHECK_ERR
    ADD_VAR_1D(elem_ss2, NC_INT, num_side_ss2) ;
    ADD_VAR_1D(side_ss2, NC_INT, num_side_ss2) ;
    ADD_VAR_1D(elem_ss5, NC_INT, num_side_ss5) ;
    ADD_VAR_1D(side_ss5, NC_INT, num_side_ss5) ;
    ADD_VAR_1D(elem_ss9, NC_INT, num_side_ss9) ;
    ADD_VAR_1D(side_ss9, NC_INT, num_side_ss9) ;
    ADD_VAR_1D(elem_ss12, NC_INT, num_side_ss12) ;
    ADD_VAR_1D(side_ss12, NC_INT, num_side_ss12) ;
    ADD_VAR_1D(elem_ss14, NC_INT, num_side_ss14) ;
    ADD_VAR_1D(side_ss14, NC_INT, num_side_ss14) ;
    ADD_VAR_1D(elem_ss15, NC_INT, num_side_ss15) ;
    ADD_VAR_1D(side_ss15, NC_INT, num_side_ss15) ;
    ADD_VAR_1D(elem_ss17, NC_INT, num_side_ss17) ;
    ADD_VAR_1D(side_ss17, NC_INT, num_side_ss17) ;
    ADD_VAR_1D(elem_ss21, NC_INT, num_side_ss21) ;
    ADD_VAR_1D(side_ss21, NC_INT, num_side_ss21) ;
    ADD_VAR_1D(elem_ss25, NC_INT, num_side_ss25) ;
    ADD_VAR_1D(side_ss25, NC_INT, num_side_ss25) ;
    ADD_VAR_1D(em_prop1, NC_INT, num_elem_maps) ;
    err = PUT_ATTR_TEXT(ncid, em_prop1_v, "name", 2, "ID"); CHECK_ERR
    ADD_VAR_2D(emap_names, NC_CHAR, num_elem_maps, len_name) ;
    err = PUT_ATTR_TEXT(ncid, emap_names_v, "_FillValue", 1, "\0"); CHECK_ERR
    ADD_VAR_1D(elem_map1, NC_INT, num_elem) ;

    err = ENDDEF(ncid); CHECK_ERR

    /* write variable eb_names */
    buf = (char*) calloc(num_el_blk * LEN_NAME, 1);
    ptr = buf;
    strcpy(ptr, "block_1");        ptr += LEN_NAME;
    strcpy(ptr, "block_1_air");    ptr += LEN_NAME;
    strcpy(ptr, "block_1_fluid");
    err = PUT_TEXT(ncid, eb_names_v, buf); CHECK_ERR
    free(buf);

    /* add a global attribute to grow header */
    err = REDEF(ncid); CHECK_ERR

    char *attr = (char*) calloc(512, 1);
    err = PUT_ATTR_TEXT(ncid, NC_GLOBAL, "attr", 512, attr);
    CHECK_ERR
    free(attr);

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
    err = PUT_TEXT(ncid, ss_names_v, buf); CHECK_ERR
    free(buf);

    /* add aother global attribute to grow header */
    err = REDEF(ncid); CHECK_ERR

    char *attr2 = (char*) calloc(512, 1);
    err = PUT_ATTR_TEXT(ncid, NC_GLOBAL, "attr2", 512, attr2);
    CHECK_ERR
    free(attr2);

    err = ENDDEF(ncid); CHECK_ERR

    /* close file */
    err = CLOSE(ncid); CHECK_ERR

    MPI_Finalize();
    return (nerrs > 0);
}

