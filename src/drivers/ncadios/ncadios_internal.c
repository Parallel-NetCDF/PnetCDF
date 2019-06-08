/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>

#include "ncadios_internal.h"

int ncadiosi_inq_varid(NC_ad* ncadp, char* name, int *id) {
    int tmp;

    if (id != NULL){
        tmp = ncadiosi_var_list_find(&(ncadp->vars), name);
        if (tmp < 0){
            DEBUG_RETURN_ERROR(NC_ENOTVAR)
        }
        *id = tmp;
    }

    return NC_NOERR;
}

int ncadiosi_inq_dimid(NC_ad* ncadp, char* name, int *id) {
    int tmp;

    if (id != NULL){
        tmp = ncadiosi_dim_list_find(&(ncadp->dims), name);
        if (tmp < 0){
            DEBUG_RETURN_ERROR(NC_EBADDIM)
        }
        *id = tmp;
    }

    return NC_NOERR;
}

int ncadiosi_def_var(NC_ad* ncadp, char* name, nc_type type, int ndim,
                        int *dimids, int *id) {
    NC_ad_var var;

    if (CHECK_NAME(name)){
        var.type = type;
        var.ndim = ndim;
        var.dimids = NCI_Malloc(SIZEOF_INT * ndim);
        memcpy(var.dimids, dimids, SIZEOF_INT * ndim);
        var.name = NCI_Malloc(strlen(name) + 1);
        strcpy(var.name, name);
        ncadiosi_att_list_init(&(var.atts));
        *id = ncadiosi_var_list_add(&(ncadp->vars), var);
    }

    return NC_NOERR;
}

int ncadiosi_def_dim(NC_ad* ncadp, char* name, int len, int *id) {
    NC_ad_dim dim;

    if (len == NC_UNLIMITED){
        *id = INT_MAX;
        return NC_NOERR;
    }

    if (CHECK_NAME(name)){
        dim.len = len;
        dim.name = NCI_Malloc(strlen(name) + 1);
        strcpy(dim.name, name);

        *id = ncadiosi_dim_list_add(&(ncadp->dims), dim);
    }

    return NC_NOERR;
}

int ncadiosi_parse_attrs(NC_ad* ncadp) {
    int i, j;
    int vid;
    char path[1024];
    char *vname = NULL;

    for (i = 0; i < ncadp->fp->nattrs; i++) {
        strcpy(path, ncadp->fp->attr_namelist[i]);

        for(j = 0; path[j] != '\0'; j++){
            if (j > 0 && path[j] == '/'){
                path[j] = '\0';
                if (path[0] == '/'){
                    vname = path + 1;
                }
                else{
                    vname = path;
                }

                vid = ncadiosi_var_list_find(&(ncadp->vars), vname);
                if (vid > -1){
                    ncadiosi_att_list_add(&(ncadp->vars.data[vid].atts), i);
                }

                break;
            }
        }
    }

    return NC_NOERR;
}

int ncadiosi_parse_rec_dim(NC_ad *ncadp) {
    int err;
    int i;
    char name[128];

    for(i = 0; i < ncadp->vars.cnt; i++){
        if (ncadp->vars.data[i].dimids[0] == INT_MAX){
            ADIOS_VARINFO * v;

            v = adios_inq_var_byid (ncadp->fp, i);
            if (v == NULL){
                err = ncmpii_error_adios2nc(adios_errno, "inq_var");
                DEBUG_RETURN_ERROR(err);
            }

            sprintf(name, "var_%d_timesteps", i);
            err = ncadiosi_def_dim(ncadp, name, v->nsteps,
                                    ncadp->vars.data[i].dimids);
            if (err != NC_NOERR){
                DEBUG_RETURN_ERROR(err)
            }
        }
    }

    return NC_NOERR;
}

int ncadiosi_parse_header_readall (NC_ad *ncadp) {
    int err;
    int i, j;
    int varid;
    int *dimids = NULL;
    int recdimid = -1;
    int maxndim = 0;
    char name[1024];

    /* For all variables */
    for (i = 0; i < ncadp->fp->nvars; i++) {
        ADIOS_VARINFO * v;

        v = adios_inq_var_byid (ncadp->fp, i);
        if (v == NULL){
            err = ncmpii_error_adios2nc(adios_errno, "inq_var");
            DEBUG_RETURN_ERROR(err);
        }

        if (maxndim < v->ndim + 1){
            maxndim = v->ndim + 1;
            if (dimids != NULL){
                NCI_Free(dimids);
            }
            dimids = (int*)NCI_Malloc(SIZEOF_INT * maxndim);
        }

        /* Record every dimensions */
        if (v->nsteps > 1){
            if (recdimid < 0){
                err = ncadiosi_def_dim(ncadp, name, NC_UNLIMITED, &recdimid);
                if (err != NC_NOERR){
                    DEBUG_RETURN_ERROR(err)
                }
            }
            dimids[0] = recdimid;
        }

        for (j = 1; j <= v->ndim; j++){
            sprintf(name, "var_%d_dim_%d", i, j);
            err = ncadiosi_def_dim(ncadp, name, v->dims[j - 1], dimids + j);
            if (err != NC_NOERR){
                DEBUG_RETURN_ERROR(err)
            }
        }

        /* Record variable */
        if (v->nsteps > 1){
            err = ncadiosi_def_var(ncadp, ncadp->fp->var_namelist[i],
                                    ncadios_to_nc_type(v->type), v->ndim + 1,
                                    dimids, &varid);
        }
        else{
            err = ncadiosi_def_var(ncadp, ncadp->fp->var_namelist[i],
                                    ncadios_to_nc_type(v->type), v->ndim,
                                    dimids + 1, &varid);
        }
        if (err != NC_NOERR){
            DEBUG_RETURN_ERROR(err)
        }

        adios_free_varinfo (v);
    } /* variables */

    return NC_NOERR;
}
