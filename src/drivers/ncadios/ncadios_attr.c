/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_inq_attname() : dispatcher->inq_attname()
 * ncmpi_inq_attid()   : dispatcher->inq_attid()
 * ncmpi_inq_att()     : dispatcher->inq_att()
 * ncmpi_rename_att()  : dispatcher->inq_rename_att()
 * ncmpi_copy_att()    : dispatcher->inq_copy_att()
 * ncmpi_del_att()     : dispatcher->inq_del_att()
 * ncmpi_get_att()     : dispatcher->inq_get_att()
 * ncmpi_put_att()     : dispatcher->inq_put_arr()
 *
 */

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


int
ncadios_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Only global attr is defined */
    /*
    if (varid >= 0){
        DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
    }
    */

    if (name != NULL){
        if (varid == NC_GLOBAL){
            if (attid >= ncadp->atts.cnt){
                DEBUG_RETURN_ERROR(NC_EINVAL);
            }
            strcpy(name, ncadp->atts.data[attid].name);
        }
        else{
            if (varid >= ncadp->vars.cnt){
                DEBUG_RETURN_ERROR(NC_EINVAL);
            }
            if (attid >= ncadp->vars.data[varid].atts.cnt){
                DEBUG_RETURN_ERROR(NC_EINVAL);
            }
            strcpy(name, ncadp->vars.data[varid].atts.data[attid].name);
        }
    }

    return NC_NOERR;
}

int
ncadios_inq_attid(void       *ncdp,
                int         varid,
                const char *name,
                int        *attidp)
{
    int err;
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (varid == NC_GLOBAL){
        for(i = 0; i < ncadp->fp->nattrs; i++){
            if (strcmp(name, ncadp->fp->attr_namelist[i]) == 0){
                if (attidp != NULL){
                    *attidp = i;
                }
                break;
            }
        }

        // Name not found
        if (i >= ncadp->fp->nattrs){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
    }
    else{
        char attname[1024];
        NC_ad_var var;

        if (varid > ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        var = ncadp->vars.data[varid];

        sprintf(attname, "/%s/%s", var.name, name);


        for (i = 0; i < ncadp->fp->nattrs; i++){
            if (strcmp(ncadp->fp->attr_namelist[i], attname) == 0){
                if (attidp != NULL){
                    *attidp = ncadiosi_att_list_find(&(var.atts), i);
                }
                break;
            }
        }
        
        // Name not found
        if (i >= ncadp->fp->nattrs){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
    }

    return NC_NOERR;

    //return ncadiosi_inq_attid(ncadp, varid, name, attidp);
}

int
ncadios_inq_att(void       *ncdp,
              int         varid,
              const char *name,
              nc_type    *datatypep,
              MPI_Offset *lenp)
{
    int err;
    int attid;
    NC_ad *ncadp = (NC_ad*)ncdp;
    NC_ad_att att;
    
    if (varid == NC_GLOBAL){
        err = adios_get_attr(ncadp->fp, name, &atype, &asize, &adata);
    }
    else{
        char attname[1024];
        if (varid >= ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        sprintf(attname, "/%s/%s", ncadp->vars.data[varid].name, name);
        err = adios_get_attr(ncadp->fp, attname, &atype, &asize, &adata);
    }

    if (err != 0){
        //TODO: translate adios error
        return err;
    }

    tsize = adios_type_size(atype, adata);

    if (datatypep != NULL){
        *datatypep = ncadios_to_nc_type(atype);
    }

    if (lenp != NULL){
        if (atype == adios_string){
            *lenp = (MPI_Offset)asize;
        }
        else{
            *lenp = (MPI_Offset)asize / tsize;
        }
    }

    return NC_NOERR;
}

int
ncadios_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}


int
ncadios_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    int err;
    NC_ad *ncadp_in  = (NC_ad*)ncdp_in;
    NC_ad *ncadp_out = (NC_ad*)ncdp_out;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadios_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    int err;
    int attid, esize;
    NC_ad *ncadp = (NC_ad*)ncdp;
    NC_ad_att att;
    MPI_Datatype xtype;
    
    if (varid == NC_GLOBAL){
        err = adios_get_attr(ncadp->fp, name, &atype, &asize, &adata);
    }
    else{
        char attname[1024];
        if (varid >= ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        sprintf(attname, "/%s/%s", ncadp->vars.data[varid].name, name);
        err = adios_get_attr(ncadp->fp, attname, &atype, &asize, &adata);
    }
    if (err != 0){
        //TODO: translate adios error
        return err;
    }

    xtype = ncadios_to_mpi_type(atype);
    if (xtype != itype){
        esize = adios_type_size (atype, adata);
        nelems = asize / esize;
        ncadiosiconvert(adata, buf, xtype, itype, nelems);
    }
    else{
        memcpy(buf, adata, asize);
    }

    /*
    if (varid == NC_GLOBAL){
        attid = ncadiosi_att_list_find(&(ncadp->atts), name);
        if (attid < 0){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        att = ncadp->atts.data[attid];
    }
    else{
        if (varid >= ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        attid = ncadiosi_att_list_find(&(ncadp->vars.data[varid].atts), name);
        if (attid < 0){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        att = ncadp->vars.data[varid].atts.data[attid];
    }

    xtype = ncadios_nc_to_mpi_type(att.type);

    if (xtype != itype){
        ncadiosiconvert(att.data, buf, xtype, itype, att.len);
    }
    else{
        MPI_Type_size(xtype, &esize);
        memcpy(buf, att.data, att.len * esize);
    }

    return NC_NOERR;
}

int
ncadios_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void   *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}
