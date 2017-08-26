/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* This program tests whether the correct error codes can be returned when
 * using NULL arguments for start, count, stride, or imap
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <netcdf.h>

#define EXP_ERR(exp,msg) { \
    if (err != exp) { \
        nerrs++; \
        fprintf(stderr, "Error at line %d in %s: %s expect (%s) but got (%s)\n", \
                __LINE__, __FILE__, msg, \
                nc_strerror(exp), nc_strerror(err)); \
    } \
}

int main(int argc, char **argv)
{
    char filename[256];
    int err, nerrs=0, ncid, dimid[2], varid;
    int buf[100];
    size_t start[2], count[2];
    ptrdiff_t stride[2], imap[2];

    if (argc > 2) {
        printf("Usage: %s [filename]\n",argv[0]);
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    err = nc_create(filename, NC_CLOBBER, &ncid);
    EXP_ERR(NC_NOERR, "create")
    err = nc_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]);
    EXP_ERR(NC_NOERR,"def_dim Y")
    err = nc_def_dim(ncid, "X", 10, &dimid[1]);
    EXP_ERR(NC_NOERR,"def_dim X")
    err = nc_def_var(ncid, "var", NC_INT, 2, dimid, &varid);
    EXP_ERR(NC_NOERR,"def_var")
    err = nc_enddef(ncid);
    EXP_ERR(NC_NOERR,"enddef")

     start[0] =  start[1] = 0;
     count[0] =  count[1] = 1;
    stride[0] = stride[1] = 1;
      imap[0] =   imap[1] = 1;

    /*---- test put_var1 ---- */
    err = nc_put_var1_int(ncid, varid, start, buf);
    EXP_ERR(NC_NOERR, "put_var1")

    err = nc_put_var1_int(ncid, varid, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_var1 start=NULL")

    /*---- test put_vara ---- */
    err = nc_put_vara_int(ncid, varid, start, count, buf);
    EXP_ERR(NC_NOERR, "put_vara")

    err = nc_put_vara_int(ncid, varid, NULL, count, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vara start=NULL")

    err = nc_put_vara_int(ncid, varid, start, NULL, buf);
    EXP_ERR(NC_EEDGE, "put_vara count=NULL")

    err = nc_put_vara_int(ncid, varid, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vara start=count=NULL")

    /*---- test put_vars ---- */
    err = nc_put_vars_int(ncid, varid, start, count, stride, buf);
    EXP_ERR(NC_NOERR, "put_vars")

    err = nc_put_vars_int(ncid, varid, NULL, count, stride, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vars start=NULL")

    err = nc_put_vars_int(ncid, varid, start, NULL, stride, buf);
    EXP_ERR(NC_EEDGE, "put_vars count=NULL")

    err = nc_put_vars_int(ncid, varid, start, count, NULL, buf);
    EXP_ERR(NC_NOERR, "put_vars stride=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, NULL, stride, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vars start=count=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, count, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vars start=stride=NULL")

    err = nc_put_vars_int(ncid, varid, start, NULL, NULL, buf);
    EXP_ERR(NC_EEDGE, "put_vars count=stride=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_vars start=count=stride=NULL")

    /*---- test put_varm ---- */
    err = nc_put_varm_int(ncid, varid, start, count, stride, imap, buf);
    EXP_ERR(NC_NOERR, "put_varm")

    err = nc_put_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    EXP_ERR(NC_EEDGE, "put_varm count=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, NULL, imap, buf);
    EXP_ERR(NC_NOERR, "put_varm stride=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, stride, NULL, buf);
    EXP_ERR(NC_NOERR, "put_varm imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=count=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=stride=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    EXP_ERR(NC_EEDGE, "put_varm count=stride=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    EXP_ERR(NC_EEDGE, "put_varm count=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    EXP_ERR(NC_NOERR, "put_varm stride=imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=count=stride=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=count=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EEDGE, "put_varm count=stride=imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "put_varm start=count=stride=imap=NULL")

    /*---- test get_var1 ---- */
    err = nc_get_var1_int(ncid, varid, start, buf);
    EXP_ERR(NC_NOERR, "get_var1")

    err = nc_get_var1_int(ncid, varid, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_var1 start=NULL")

    /*---- test get_vara ---- */
    err = nc_get_vara_int(ncid, varid, start, count, buf);
    EXP_ERR(NC_NOERR, "get_vara")

    err = nc_get_vara_int(ncid, varid, NULL, count, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vara start=NULL")

    err = nc_get_vara_int(ncid, varid, start, NULL, buf);
    EXP_ERR(NC_EEDGE, "get_vara count=NULL")

    err = nc_get_vara_int(ncid, varid, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vara start=count=NULL")

    /*---- test get_vars ---- */
    err = nc_get_vars_int(ncid, varid, start, count, stride, buf);
    EXP_ERR(NC_NOERR, "get_vars")

    err = nc_get_vars_int(ncid, varid, NULL, count, stride, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vars start=NULL")

    err = nc_get_vars_int(ncid, varid, start, NULL, stride, buf);
    EXP_ERR(NC_EEDGE, "get_vars count=NULL")

    err = nc_get_vars_int(ncid, varid, start, count, NULL, buf);
    EXP_ERR(NC_NOERR, "get_vars stride=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, NULL, stride, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vars start=count=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, count, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vars start=stride=NULL")

    err = nc_get_vars_int(ncid, varid, start, NULL, NULL, buf);
    EXP_ERR(NC_EEDGE, "get_vars count=stride=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_vars start=count=stride=NULL")

    /*---- test get_varm ---- */
    err = nc_get_varm_int(ncid, varid, start, count, stride, imap, buf);
    EXP_ERR(NC_NOERR, "get_varm")

    err = nc_get_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    EXP_ERR(NC_EEDGE, "get_varm count=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, NULL, imap, buf);
    EXP_ERR(NC_NOERR, "get_varm stride=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, stride, NULL, buf);
    EXP_ERR(NC_NOERR, "get_varm imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=count=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=stride=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    EXP_ERR(NC_EEDGE, "get_varm count=stride=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    EXP_ERR(NC_EEDGE, "get_varm count=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    EXP_ERR(NC_NOERR, "get_varm stride=imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=count=stride=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=count=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EEDGE, "get_varm count=stride=imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    EXP_ERR(NC_EINVALCOORDS, "get_varm start=count=stride=imap=NULL")

    err = nc_close(ncid);
    EXP_ERR(NC_NOERR, "close")

    return (nerrs > 0);
}

