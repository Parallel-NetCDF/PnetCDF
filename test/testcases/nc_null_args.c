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

static char* nc_err_code_name(int err);

#define ERR(exp,msg) { \
    if (err != exp) { \
        nerrs++; \
        fprintf(stderr, "Error at line %d: %s expect %s but got %s\n", \
                __LINE__, msg, \
                nc_err_code_name(exp), nc_err_code_name(err)); \
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
        return 0;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    err = nc_create(filename, NC_CLOBBER, &ncid);
    ERR(NC_NOERR, "create")
    err = nc_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]);
    ERR(NC_NOERR,"def_dim Y")
    err = nc_def_dim(ncid, "X", 10, &dimid[1]);
    ERR(NC_NOERR,"def_dim X")
    err = nc_def_var(ncid, "var", NC_INT, 2, dimid, &varid);
    ERR(NC_NOERR,"def_var")
    err = nc_enddef(ncid);
    ERR(NC_NOERR,"enddef")

     start[0] =  start[1] = 0;
     count[0] =  count[1] = 1;
    stride[0] = stride[1] = 1;
      imap[0] =   imap[1] = 1;

    /*---- test put_var1 ---- */
    err = nc_put_var1_int(ncid, varid, start, buf);
    ERR(NC_NOERR, "put_var1")

    err = nc_put_var1_int(ncid, varid, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_var1 start=NULL")

    /*---- test put_vara ---- */
    err = nc_put_vara_int(ncid, varid, start, count, buf);
    ERR(NC_NOERR, "put_vara")

    err = nc_put_vara_int(ncid, varid, NULL, count, buf);
    ERR(NC_EINVALCOORDS, "put_vara start=NULL")

    err = nc_put_vara_int(ncid, varid, start, NULL, buf);
    ERR(NC_EEDGE, "put_vara count=NULL")

    err = nc_put_vara_int(ncid, varid, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_vara start=count=NULL")

    /*---- test put_vars ---- */
    err = nc_put_vars_int(ncid, varid, start, count, stride, buf);
    ERR(NC_NOERR, "put_vars")

    err = nc_put_vars_int(ncid, varid, NULL, count, stride, buf);
    ERR(NC_EINVALCOORDS, "put_vars start=NULL")

    err = nc_put_vars_int(ncid, varid, start, NULL, stride, buf);
    ERR(NC_EEDGE, "put_vars count=NULL")

    err = nc_put_vars_int(ncid, varid, start, count, NULL, buf);
    ERR(NC_NOERR, "put_vars stride=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, NULL, stride, buf);
    ERR(NC_EINVALCOORDS, "put_vars start=count=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, count, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_vars start=stride=NULL")

    err = nc_put_vars_int(ncid, varid, start, NULL, NULL, buf);
    ERR(NC_EEDGE, "put_vars count=stride=NULL")

    err = nc_put_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_vars start=count=stride=NULL")

    /*---- test put_varm ---- */
    err = nc_put_varm_int(ncid, varid, start, count, stride, imap, buf);
    ERR(NC_NOERR, "put_varm")

    err = nc_put_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    ERR(NC_EEDGE, "put_varm count=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, NULL, imap, buf);
    ERR(NC_NOERR, "put_varm stride=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, stride, NULL, buf);
    ERR(NC_NOERR, "put_varm imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=count=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=stride=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    ERR(NC_EEDGE, "put_varm count=stride=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    ERR(NC_EEDGE, "put_varm count=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    ERR(NC_NOERR, "put_varm stride=imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=count=stride=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=count=imap=NULL")

    err = nc_put_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    ERR(NC_EEDGE, "put_varm count=stride=imap=NULL")

    err = nc_put_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "put_varm start=count=stride=imap=NULL")

    /*---- test get_var1 ---- */
    err = nc_get_var1_int(ncid, varid, start, buf);
    ERR(NC_NOERR, "get_var1")

    err = nc_get_var1_int(ncid, varid, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_var1 start=NULL")

    /*---- test get_vara ---- */
    err = nc_get_vara_int(ncid, varid, start, count, buf);
    ERR(NC_NOERR, "get_vara")

    err = nc_get_vara_int(ncid, varid, NULL, count, buf);
    ERR(NC_EINVALCOORDS, "get_vara start=NULL")

    err = nc_get_vara_int(ncid, varid, start, NULL, buf);
    ERR(NC_EEDGE, "get_vara count=NULL")

    err = nc_get_vara_int(ncid, varid, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_vara start=count=NULL")

    /*---- test get_vars ---- */
    err = nc_get_vars_int(ncid, varid, start, count, stride, buf);
    ERR(NC_NOERR, "get_vars")

    err = nc_get_vars_int(ncid, varid, NULL, count, stride, buf);
    ERR(NC_EINVALCOORDS, "get_vars start=NULL")

    err = nc_get_vars_int(ncid, varid, start, NULL, stride, buf);
    ERR(NC_EEDGE, "get_vars count=NULL")

    err = nc_get_vars_int(ncid, varid, start, count, NULL, buf);
    ERR(NC_NOERR, "get_vars stride=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, NULL, stride, buf);
    ERR(NC_EINVALCOORDS, "get_vars start=count=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, count, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_vars start=stride=NULL")

    err = nc_get_vars_int(ncid, varid, start, NULL, NULL, buf);
    ERR(NC_EEDGE, "get_vars count=stride=NULL")

    err = nc_get_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_vars start=count=stride=NULL")

    /*---- test get_varm ---- */
    err = nc_get_varm_int(ncid, varid, start, count, stride, imap, buf);
    ERR(NC_NOERR, "get_varm")

    err = nc_get_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    ERR(NC_EEDGE, "get_varm count=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, NULL, imap, buf);
    ERR(NC_NOERR, "get_varm stride=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, stride, NULL, buf);
    ERR(NC_NOERR, "get_varm imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=count=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=stride=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    ERR(NC_EEDGE, "get_varm count=stride=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    ERR(NC_EEDGE, "get_varm count=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    ERR(NC_NOERR, "get_varm stride=imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=count=stride=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=count=imap=NULL")

    err = nc_get_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    ERR(NC_EEDGE, "get_varm count=stride=imap=NULL")

    err = nc_get_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    ERR(NC_EINVALCOORDS, "get_varm start=count=stride=imap=NULL")

    err = nc_close(ncid);
    ERR(NC_NOERR, "close")

    return 0;
}

static
char* nc_err_code_name(int err)
{
    static char unknown_str[32];

    if (err > 0) { /* system error */
        const char *cp = (const char *) strerror(err);
        if (cp == NULL)
            sprintf(unknown_str,"Unknown error code %d",err);
        else
            sprintf(unknown_str,"Error code %d (%s)",err,cp);
        return unknown_str;
    }

    switch (err) {
        case (NC_NOERR):			return "NC_NOERR";
        case (NC_EBADID):			return "NC_EBADID";
        case (NC_ENFILE):			return "NC_ENFILE";
        case (NC_EEXIST):			return "NC_EEXIST";
        case (NC_EINVAL):			return "NC_EINVAL";
        case (NC_EPERM):			return "NC_EPERM";
        case (NC_ENOTINDEFINE):			return "NC_ENOTINDEFINE";
        case (NC_EINDEFINE):			return "NC_EINDEFINE";
        case (NC_EINVALCOORDS):			return "NC_EINVALCOORDS";
        case (NC_EMAXDIMS):			return "NC_EMAXDIMS";
        case (NC_ENAMEINUSE):			return "NC_ENAMEINUSE";
        case (NC_ENOTATT):			return "NC_ENOTATT";
        case (NC_EMAXATTS):			return "NC_EMAXATTS";
        case (NC_EBADTYPE):			return "NC_EBADTYPE";
        case (NC_EBADDIM):			return "NC_EBADDIM";
        case (NC_EUNLIMPOS):			return "NC_EUNLIMPOS";
        case (NC_EMAXVARS):			return "NC_EMAXVARS";
        case (NC_ENOTVAR):			return "NC_ENOTVAR";
        case (NC_EGLOBAL):			return "NC_EGLOBAL";
        case (NC_ENOTNC):			return "NC_ENOTNC";
        case (NC_ESTS):				return "NC_ESTS";
        case (NC_EMAXNAME):			return "NC_EMAXNAME";
        case (NC_EUNLIMIT):			return "NC_EUNLIMIT";
        case (NC_ENORECVARS):			return "NC_ENORECVARS";
        case (NC_ECHAR):			return "NC_ECHAR";
        case (NC_EEDGE):			return "NC_EEDGE";
        case (NC_ESTRIDE):			return "NC_ESTRIDE";
        case (NC_EBADNAME):			return "NC_EBADNAME";
        case (NC_ERANGE):			return "NC_ERANGE";
        case (NC_ENOMEM):			return "NC_ENOMEM";
        case (NC_EVARSIZE):			return "NC_EVARSIZE";
        case (NC_EDIMSIZE):			return "NC_EDIMSIZE";
        case (NC_ETRUNC):			return "NC_ETRUNC";
        case (NC_EAXISTYPE):			return "NC_EAXISTYPE";
        case (NC_EDAP):				return "NC_EDAP";
        case (NC_ECURL):			return "NC_ECURL";
        case (NC_EIO):				return "NC_EIO";
        case (NC_ENODATA):			return "NC_ENODATA";
        case (NC_EDAPSVC):			return "NC_EDAPSVC";
        case (NC_EDAS):				return "NC_EDAS";
        case (NC_EDDS):				return "NC_EDDS";
        case (NC_EDATADDS):			return "NC_EDATADDS";
        case (NC_EDAPURL):			return "NC_EDAPURL";
        case (NC_EDAPCONSTRAINT):		return "NC_EDAPCONSTRAINT";
        case (NC_ETRANSLATION):			return "NC_ETRANSLATION";
        case (NC_EACCESS):			return "NC_EACCESS";
        case (NC_EAUTH):			return "NC_EAUTH";
        case (NC_ENOTFOUND):			return "NC_ENOTFOUND";
        case (NC_ECANTREMOVE):			return "NC_ECANTREMOVE";
        case (NC_EHDFERR):			return "NC_EHDFERR";
        case (NC_ECANTREAD):			return "NC_ECANTREAD";
        case (NC_ECANTWRITE):			return "NC_ECANTWRITE";
        case (NC_ECANTCREATE):			return "NC_ECANTCREATE";
        case (NC_EFILEMETA):			return "NC_EFILEMETA";
        case (NC_EDIMMETA):			return "NC_EDIMMETA";
        case (NC_EATTMETA):			return "NC_EATTMETA";
        case (NC_EVARMETA):			return "NC_EVARMETA";
        case (NC_ENOCOMPOUND):			return "NC_ENOCOMPOUND";
        case (NC_EATTEXISTS):			return "NC_EATTEXISTS";
        case (NC_ENOTNC4):			return "NC_ENOTNC4";
        case (NC_ESTRICTNC3):			return "NC_ESTRICTNC3";
        case (NC_ENOTNC3):			return "NC_ENOTNC3";
        case (NC_ENOPAR):			return "NC_ENOPAR";
        case (NC_EPARINIT):			return "NC_EPARINIT";
        case (NC_EBADGRPID):			return "NC_EBADGRPID";
        case (NC_EBADTYPID):			return "NC_EBADTYPID";
        case (NC_ETYPDEFINED):			return "NC_ETYPDEFINED";
        case (NC_EBADFIELD):			return "NC_EBADFIELD";
        case (NC_EBADCLASS):			return "NC_EBADCLASS";
        case (NC_EMAPTYPE):			return "NC_EMAPTYPE";
        case (NC_ELATEFILL):			return "NC_ELATEFILL";
        case (NC_ELATEDEF):			return "NC_ELATEDEF";
        case (NC_EDIMSCALE):			return "NC_EDIMSCALE";
        case (NC_ENOGRP):			return "NC_ENOGRP";
        case (NC_ESTORAGE):			return "NC_ESTORAGE";
        case (NC_EBADCHUNK):			return "NC_EBADCHUNK";
        case (NC_ENOTBUILT):			return "NC_ENOTBUILT";
        case (NC_EDISKLESS):			return "NC_EDISKLESS";
        case (NC_ECANTEXTEND):			return "NC_ECANTEXTEND";
        case (NC_EMPI):				return "NC_EMPI";
        // case (NC_EURL):				return "NC_EURL";
        // case (NC_ECONSTRAINT):			return "NC_ECONSTRAINT";
        default:
              sprintf(unknown_str,"Unknown code %d",err);
    }
    return unknown_str;
}

