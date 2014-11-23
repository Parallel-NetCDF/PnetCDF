/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can detect file header inconsistency and
 * overwrite the inconsistent header with root's
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>

#define ERR_EXP(e, exp) {if (e != exp) { printf("Error (line %d): expecting error code %d but got %d\n", __LINE__, exp, e); nerr++; }}
#define ERR_EXP2(e, exp1, exp2) {if (e != exp1 && e != exp2) { printf("Error (line %d): expecting error code %d or %d but got %d\n", __LINE__, exp1, exp2, e); nerr++; }}

#define ERR {if(err!=NC_NOERR) {printf("Error(%d) at line %d: %s\n",err,__LINE__,ncmpi_strerror(err)); nerr++; }}

/*----< test_open_mode() >----------------------------------------------------*/
static
int test_open_mode(char *filename, int safe_mode)
{
    int err, rank, ncid, cmode, omode, nerr=0;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);

    /* Test inconsistent cmode -----------------------------------------------*/
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;
    if (rank == 0) cmode = NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid);
    /* if safe_mode is on, then we expect all non-root ranks to print a warning
     * message "inconsistent file create mode, overwrite with root's" and error
     * code NC_EMULTIDEFINE_OMODE from non-root processes.
     * if safe_mode is off, then no error code should be returned.
     */
    if (safe_mode && rank > 0) ERR_EXP(err, NC_EMULTIDEFINE_OMODE)
    else                       ERR

    err = ncmpi_close(ncid);
    /* if safe_mode is on, then no error code should be returned.
     * if safe_mode is off, then we expect error code NC_EMULTIDEFINE_OMODE
     * from all non-root processes and NC_EMULTIDEFINE from root process.
     */
    if (safe_mode) ERR
    else           ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_OMODE)

    int format;
    err = ncmpi_inq_file_format(filename, &format); ERR
    if (format != 1) {
        printf("Error (line %d): output file should be in CDF-1 format\n",__LINE__);
        nerr++;
    }

    /* Test inconsistent omode -----------------------------------------------*/
    omode = NC_WRITE;
    if (rank == 0) omode = NC_NOWRITE;
    err = ncmpi_open(comm, filename, omode, info, &ncid);
    if (safe_mode) {
        /* in safe_mode, we expect all non-root ranks to print a warning message
         * "inconsistent file open mode, overwrite with root's" and error code
         * code NC_EMULTIDEFINE_OMODE and no error on root process.
         */
        /* expected errors: NC_EMULTIDEFINE or NC_EMULTIDEFINE_OMODE */
        if (rank > 0) ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_OMODE)

        /* in safe mode, open mode inconsistent is not a fatal error, file is
         * still opened, with root's omode overwriting others. Hence, once the
         * test is done, we need to close the file */
        err = ncmpi_close(ncid); ERR
    }
    else {
        /* expected errors: NC_EMULTIDEFINE_OMODE or NC_EMULTIDEFINE_FNC_ARGS */
        ERR_EXP2(err, NC_EMULTIDEFINE_OMODE, NC_EMULTIDEFINE_FNC_ARGS)

        /* When not in safe mode, the inconsistent omode will be passed to
         * MPI_File_open(). MPI-IO should return error class MPI_ERR_NOT_SAME
         * which will be translated to NC_EMULTIDEFINE_FNC_ARGS in PnetCDF.
         * If MPI-IO reports error class MPI_ERR_AMODE instead, then it will be
         * translated to NC_EMULTIDEFINE_OMODE in PnetCDF. In any case, the file
         * will not be opened by MPI-IO and hence we need not close the file.
         */
    }
    return nerr;
}

/*----< test_dim() >----------------------------------------------------------*/
static
int test_dim(char *filename)
{
    int err, rank, ncid, cmode, ndims, dimid1, dimid2, dimid3, nerr=0;
    MPI_Offset dimlen;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistency on number of dimensions ----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "x", 100, &dimid1); ERR
    if (rank == 0) {
        err = ncmpi_def_dim(ncid, "y", 100, &dimid2); ERR
    }
    err = ncmpi_enddef(ncid);
    /* expected errors: NC_EMULTIDEFINE or NC_EMULTIDEFINE_DIM_NUM */
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_DIM_NUM)

    err = ncmpi_inq_ndims(ncid, &ndims); ERR
    if (ndims != 2) {
        printf("Error (line %d): number of dimesnions (%d) defined should be 2\n",__LINE__,ndims);
        nerr++;
    }
    /* all processes should be able to see dim "y" */
    err = ncmpi_inq_dimid(ncid, "y", &dimid2);
    if (err != NC_NOERR) {
        printf("Error (line %d): all processes should be able to see dim \"y\"\n",__LINE__);
        nerr++;
        ERR
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistency on number of dimensions ----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "x", 100, &dimid1); ERR
    if (rank > 0) {
        err = ncmpi_def_dim(ncid, "y", 100, &dimid2); ERR
    }
    err = ncmpi_enddef(ncid);
    /* expected errors: NC_EMULTIDEFINE or NC_EMULTIDEFINE_DIM_NUM */
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_DIM_NUM)

    err = ncmpi_inq_ndims(ncid, &ndims); ERR
    if (ndims != 1) {
        printf("Error (line %d): number of dimesnions (%d) defined should be 1\n",__LINE__,ndims);
        nerr++;
    }
    /* no process should be able to get dim "y" */
    err = ncmpi_inq_dimid(ncid, "y", &dimid2);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EBADDIM)

    err = ncmpi_close(ncid); ERR

    /* Test inconsistency on dimension names ---------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    if (rank == 0)
        err = ncmpi_def_dim(ncid, "y", 100, &dimid1);
    else
        err = ncmpi_def_dim(ncid, "xx", 100, &dimid1);
    ERR
    err = ncmpi_enddef(ncid);
    /* expected errors: NC_EMULTIDEFINE or NC_EMULTIDEFINE_DIM_NAME */
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_DIM_NAME)

    /* all processes should be able to get dim "y" */
    err = ncmpi_inq_dimid(ncid, "y", &dimid2); ERR

    /* no process should be able to get dim "x" */
    err = ncmpi_inq_dimid(ncid, "xx", &dimid3);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EBADDIM)

    err = ncmpi_close(ncid); ERR

    /* Test inconsistency on dimension size ----------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    if (rank == 0)
        err = ncmpi_def_dim(ncid, "x", 99, &dimid1);
    else
        err = ncmpi_def_dim(ncid, "x", 100, &dimid1);
    ERR
    err = ncmpi_enddef(ncid);
    /* expected errors: NC_EMULTIDEFINE or NC_EMULTIDEFINE_DIM_SIZE */
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_DIM_SIZE)

    /* all processes should be able to get dim "x" of size == 99 */
    err = ncmpi_inq_dimlen(ncid, dimid1, &dimlen); ERR
    if (dimlen != 99) {
        printf("Error (line %d): dimesnion size (%lld) should be 99\n",__LINE__,dimlen);
        nerr++;
    }

    err = ncmpi_close(ncid); ERR
    return nerr;
}

/*----< test_attr() >---------------------------------------------------------*/
static
int test_attr(char *filename)
{
    int err, rank, ncid, cmode, nerr=0;
    char  gattr[128];
    int   int_attr;
    float flt_attr;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistent global attribute numbers ----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    int_attr = 1;
    flt_attr = 1.0;
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gattr_1", NC_INT, 1, &int_attr);
    ERR
    if (rank == 0) {
        err = ncmpi_put_att_float(ncid, NC_GLOBAL, "gattr_2", NC_FLOAT, 1,
                                  &flt_attr); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_NUM)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute name -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    int_attr = 1;
    sprintf(gattr, "gattr_name.%d",rank);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, gattr, NC_INT, 1, &int_attr); ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_NAME)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute type -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 1, &int_attr);
    else
        err = ncmpi_put_att_float(ncid, NC_GLOBAL, "gatt", NC_FLOAT, 1, &flt_attr);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_TYPE)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    int intv[2]={1,2};
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 2, intv);
    else
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 1, intv);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_LEN)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    if (rank == 0) intv[1]=3;
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 2, intv);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_VAL)
    err = ncmpi_close(ncid); ERR

    return nerr;
}

/*----< test_var() >----------------------------------------------------------*/
static
int test_var(char *filename)
{
    int err, rank, ncid, cmode, nerr=0;
    int ndims, dimid[3], nvars, varid1, varid2, int_attr;
    float flt_attr;
    char name[128], var_attr[128];
    nc_type xtype;
    MPI_Offset dimlen;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistent variable attribute numbers --------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    int_attr = 1;
    flt_attr = 1.0;
    err = ncmpi_put_att_int(ncid, varid1, "var_attr_1", NC_INT, 1, &int_attr);
    ERR
    if (rank == 0) {
        err = ncmpi_put_att_float(ncid, varid1, "var_attr_2", NC_FLOAT, 1,
                                  &flt_attr); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_NUM)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute name -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    int_attr = 1;
    sprintf(var_attr, "var_attr_name.%d",rank);
    err = ncmpi_put_att_int(ncid, varid1, var_attr, NC_INT, 1, &int_attr); ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_NAME)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute type -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 1, &int_attr);
    else
        err = ncmpi_put_att_float(ncid, varid1, "var_att", NC_FLOAT, 1, &flt_attr);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_TYPE)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    int intv[2]={1,2};
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 2, intv);
    else
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 1, intv);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_LEN)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    if (rank == 0) intv[1]=3;
    err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 2, intv);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_ATTR_VAL)
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent number of variables ---------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); ERR
    if (rank == 0) {
        err = ncmpi_def_var(ncid, "var2", NC_INT, 1, dimid, &varid2); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_NUM)

    err = ncmpi_inq_nvars(ncid, &nvars); ERR
    if (nvars != 2) {
        printf("Error (line %d): all processes should see 2 variables\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent variable name ---------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    sprintf(name, "var.%d",rank);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, dimid, &varid1); ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_NAME)

    err = ncmpi_inq_varname(ncid, varid1, name); ERR
    if (strcmp(name, "var.0")) {
        printf("Error (line %d): all processes should see variable name: \"var.0\"\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent variable ndims --------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim0", 3, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "dim1", 2, &dimid[1]); ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 1, dimid, &varid1);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_NDIMS)

    err = ncmpi_inq_varndims(ncid, varid1, &ndims); ERR
    if (ndims != 2) {
        printf("Error (line %d): all processes should see var has 2 dimensions\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent variable type ---------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_INT, 1, dimid, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 1, dimid, &varid1);
    ERR
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_TYPE)

    err = ncmpi_inq_vartype(ncid, varid1, &xtype); ERR
    if (xtype != NC_INT) {
        printf("Error (line %d): all processes should see var is of type NC_INT\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent variable length -------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim0", 5, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "dim1", 4, &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "dim2", 3, &dimid[2]); ERR
    if (rank == 0) {
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1); ERR
    }
    else {
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid+1, &varid1); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_LEN)

    err = ncmpi_inq_vardimid(ncid, varid1, dimid); ERR
    err = ncmpi_inq_dimname(ncid, dimid[0], name); ERR
    if (strcmp(name, "dim0")) {
        printf("Error (line %d): all processes should see var's dim[0] name \"dim0\"\n",__LINE__);
        nerr++;
    }
    err = ncmpi_inq_dimname(ncid, dimid[1], name); ERR
    if (strcmp(name, "dim1")) {
        printf("Error (line %d): all processes should see var's dim[1] name \"dim1\"\n",__LINE__);
        nerr++;
    }
    err = ncmpi_inq_dimlen(ncid, dimid[0], &dimlen); ERR
    if (dimlen != 5) {
        printf("Error (line %d): all processes should see var's dim[0] len == 5\n",__LINE__);
        nerr++;
    }
    err = ncmpi_inq_dimlen(ncid, dimid[1], &dimlen); ERR
    if (dimlen != 4) {
        printf("Error (line %d): all processes should see var's dim[1] len == 4\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    /* Test inconsistent variable dimension IDs ------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "Z", 3, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y", 3, &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X", 3, &dimid[2]); ERR
    if (rank == 0) {
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid+1, &varid1); ERR
    }
    else {
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_VAR_DIMIDS)

    err = ncmpi_inq_vardimid(ncid, varid1, dimid); ERR
    err = ncmpi_inq_dimname(ncid, dimid[0], name); ERR
    if (strcmp(name, "Y")) {
        printf("Error (line %d): all processes should see var's dim[0] name \"Y\"\n",__LINE__);
        nerr++;
    }
    err = ncmpi_inq_dimname(ncid, dimid[1], name); ERR
    if (strcmp(name, "X")) {
        printf("Error (line %d): all processes should see var's dim[1] name \"X\"\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR

    return nerr;
}

/*----< test_dim_var() >------------------------------------------------------*/
static
int test_dim_var(char *filename)
{
    int i, err, rank, ncid, cmode, nerr=0;
    int ndims, dimid[3], varid1;
    char name[128], dimname[128];
    MPI_Offset dimlen;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistent of dimensions impact to variables -------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    err = ncmpi_def_dim(ncid, "dim1", 5, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "dim2", 4, &dimid[1]); ERR
    if (rank == 0) {
        err = ncmpi_def_dim(ncid, "dim3", 3, &dimid[2]); ERR
        err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimid+1, &varid1); ERR
    }
    else {
        err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimid, &varid1); ERR
    }
    err = ncmpi_enddef(ncid);
    ERR_EXP2(err, NC_EMULTIDEFINE, NC_EMULTIDEFINE_DIM_NUM)

    err = ncmpi_inq_ndims(ncid, &ndims); ERR
    if (ndims != 3) {
        printf("Error (line %d): all processes should see 3 dimensions\n",__LINE__);
        nerr++;
    }
    dimid[0] = dimid[1] = dimid[2] = -1;
    err = ncmpi_inq_vardimid(ncid, varid1, dimid); ERR
    for (i=0; i<2; i++) {
        err = ncmpi_inq_dimname(ncid, dimid[i], name); ERR
        sprintf(dimname, "dim%d", i);
        if (!strcmp(name, dimname)) {
            printf("Error (line %d): all processes should see dimid[%d] name \"%s\"\n",__LINE__,i,dimname);
            nerr++;
        }
    }
    err = ncmpi_inq_dimlen(ncid, dimid[0], &dimlen); ERR
    if (dimlen != 4) {
        printf("Error (line %d): all processes should see dimid[0] len = 4\n",__LINE__);
        nerr++;
    }
    err = ncmpi_inq_dimlen(ncid, dimid[1], &dimlen); ERR
    if (dimlen != 3) {
        printf("Error (line %d): all processes should see dimid[1] len = 3\n",__LINE__);
        nerr++;
    }
    err = ncmpi_close(ncid); ERR


    return nerr;
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *filename="testfile.nc", *mode[2] = {"0", "1"};
    int i, rank, verbose, nerr=0, sum_nerr;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    if (argc == 2) filename = argv[1];

    verbose = 0;
    for (i=verbose; i>=0; i--) {
        /* test with safe mode off and on :
         * Note even if --enable-debug is set at configure time, safe mode
         * can still be disabled by setting the environment variable
         * PNETCDF_SAFE_MODE to 0.
         */
        setenv("PNETCDF_SAFE_MODE", mode[i], 1);
        nerr += test_open_mode(filename, i);

        nerr += test_dim(filename);

        nerr += test_attr(filename);

        nerr += test_var(filename);

        nerr += test_dim_var(filename);
    }

    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Reduce(&nerr, &sum_nerr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        char cmd_str[80];
        sprintf(cmd_str, "*** TESTING C   %s for header consistency", argv[0]);
        if (sum_nerr == 0)
            printf("%-66s ------ pass\n", cmd_str);
        else
            printf("%-66s ------ failed\n", cmd_str);
    }
    MPI_Finalize();
    return 0;
}

