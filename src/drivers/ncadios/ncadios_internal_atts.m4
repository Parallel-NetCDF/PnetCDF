/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncadios_driver.h>

#include "ncadios_internal.h"

define(`PUTATT',dnl
`dnl
int ncadiosi_put_att_$2(NC_ad *ncadp, int valid, char *name, nc_type type, 
                        int len, void *value){
    double *buf;
    NC_ad_att att;
    MPI_Datatype itype;

    att.data = NCI_Malloc(len * sizeof( $2 ));
    att.len = len;
    /* Convert is not needed becuase bp2ncd never use incompetible type */
    if (ncadios_nc_to_mpi_type(type) != $1){
        printf("Warning: type mismatch\n");
    }
    memcpy(att.data, value, len * sizeof( $2 ));

    att.name = NCI_Malloc(strlen(name) + 1);
    strcpy(att.name, name);

    att.type = type;

    if (valid == NC_GLOBAL){
        ncadiosi_att_list_add(&(ncadp->atts), att);
    }
    else{
        ncadiosi_att_list_add(&(ncadp->vars.data[valid].atts), att);
    }

    return NC_NOERR;
}

')dnl
dnl
include(`foreach.m4')dnl
include(`utils.m4')dnl
dnl
foreach(`dt', (`(`MPI_CHAR', `uchar')', dnl
                `(`MPI_SIGNED_CHAR', `schar')', dnl
                `(`MPI_SHORT', `short')', dnl
                `(`MPI_INT', `int')', dnl
                `(`MPI_UNSIGNED', `long')', dnl
                `(`MPI_FLOAT', `float')', dnl
                `(`MPI_DOUBLE', `double')', dnl
                ), `PUTATT(translit(dt, `()'), $2)')dnl

int ncadiosi_put_att_text(NC_ad *ncadp, int valid, char *name, int len, 
                            void *value){
    double *buf;
    NC_ad_att att;
    MPI_Datatype itype;

    att.data = NCI_Malloc(len * SIZEOF_CHAR);
    att.len = len;
    /* Convert is not needed becuase bp2ncd never use incompetible type */
    memcpy(att.data, value, len * SIZEOF_CHAR);

    att.name = NCI_Malloc(strlen(name) + 1);
    strcpy(att.name, name);

    att.type = NC_CHAR;
    
    if (valid == NC_GLOBAL){
        ncadiosi_att_list_add(&(ncadp->atts), att);
    }
    else{
        ncadiosi_att_list_add(&(ncadp->vars.data[valid].atts), att);
    }

    return NC_NOERR;
}