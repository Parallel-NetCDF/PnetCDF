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
#include <ncadios_internal.h>


int
ncadios_inq_attname(void *ncdp,
                  int   varid,
                  int   attid,
                  char *name)
{
    NC_ad *ncadp = (NC_ad*)ncdp;

    if (varid == NC_GLOBAL){
        if (attid >= ncadp->fp->nattrs){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }

<<<<<<< HEAD
        if (name != NULL){
            strcpy(name, ncadp->fp->attr_namelist[attid]);
        }
    }
    else{
        NC_ad_var var;

        if (varid > ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        var = ncadp->vars.data[varid];

        if (attid >= var.atts.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }

        if (name != NULL){
            strcpy(name, ncadp->fp->attr_namelist[attid] + strlen(var.name) + 1);
=======
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
>>>>>>> e9cafd8... debug
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
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 23ad299... parse var attributes
    if (varid == NC_GLOBAL){
        for(i = 0; i < ncadp->fp->nattrs; i++){
            if (strcmp(name, ncadp->fp->attr_namelist[i]) == 0){
                if (attidp != NULL){
                    *attidp = i;
                }
                break;
            }
        }

<<<<<<< HEAD
        /* Name not found */
=======
        // Name not found
>>>>>>> 23ad299... parse var attributes
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

<<<<<<< HEAD
        /* ADIOS read expose all attributes as global attributes
         * Variable attributes are represented by path starting with variable 
         * name
         * These variables can be accessed as either global attributes using 
         * full path or as variable attributes using only attribute name
         */
        sprintf(attname, "%s/%s", var.name, name);
=======
        sprintf(attname, "/%s/%s", var.name, name);

>>>>>>> 23ad299... parse var attributes

        for (i = 0; i < ncadp->fp->nattrs; i++){
            if (strcmp(ncadp->fp->attr_namelist[i], attname) == 0){
                if (attidp != NULL){
                    *attidp = ncadiosi_att_list_find(&(var.atts), i);
                }
                break;
            }
        }
        
<<<<<<< HEAD
        /* Name not found */
=======
        // Name not found
>>>>>>> 23ad299... parse var attributes
        if (i >= ncadp->fp->nattrs){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
    }

    return NC_NOERR;
=======
    return ncadiosi_inq_attid(ncadp, varid, name, attidp);
>>>>>>> e9cafd8... debug
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
<<<<<<< HEAD
<<<<<<< HEAD
        err = adios_get_attr(ncadp->fp, name, &atype, &asize, &adata);
=======
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
>>>>>>> 23ad299... parse var attributes
    }
    else{
        char attname[1024];
        if (varid >= ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        
        /* ADIOS read expose all attributes as global attributes
         * Variable attributes are represented by path starting with variable 
         * name
         * These variables can be accessed as either global attributes using 
         * full path or as variable attributes using only attribute name
         */
        sprintf(attname, "/%s/%s", ncadp->vars.data[varid].name, name);
        err = adios_get_attr(ncadp->fp, attname, &atype, &asize, &adata);
    }
<<<<<<< HEAD
    if (err != 0){
        err = ncmpii_error_adios2nc(adios_errno, "get_attr");
        DEBUG_RETURN_ERROR(err);
=======
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
>>>>>>> e9cafd8... debug
    }

    if (datatypep != NULL){
        *datatypep = att.type;
    }

    if (lenp != NULL){
        *lenp = att.len;
    }

    free(adata);
=======
>>>>>>> 23ad299... parse var attributes

    return NC_NOERR;
}

int
ncadios_rename_att(void       *ncdp,
                 int         varid,
                 const char *name,
                 const char *newname)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}


int
ncadios_copy_att(void       *ncdp_in,
               int         varid_in,
               const char *name,
               void       *ncdp_out,
               int         varid_out)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}

int
ncadios_del_att(void       *ncdp,
              int         varid,
              const char *name)
{
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
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
    MPI_Datatype xtype;
    
<<<<<<< HEAD
=======
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
>>>>>>> 23ad299... parse var attributes
    if (varid == NC_GLOBAL){
        err = adios_get_attr(ncadp->fp, name, &atype, &asize, &adata);
    }
    else{
        char attname[1024];
        if (varid >= ncadp->vars.cnt){
            DEBUG_RETURN_ERROR(NC_EINVAL);
        }
        /* ADIOS read expose all attributes as global attributes
         * Variable attributes are represented by path starting with variable 
         * name
         * These variables can be accessed as either global attributes using 
         * full path or as variable attributes using only attribute name
         */
        sprintf(attname, "/%s/%s", ncadp->vars.data[varid].name, name);
        err = adios_get_attr(ncadp->fp, attname, &atype, &asize, &adata);
    }
    if (err != 0){
        err = ncmpii_error_adios2nc(adios_errno, "Open");
        DEBUG_RETURN_ERROR(err);
    }

    /* PnetCDF allow accessing attributes of different type
     * Check if we need to convert
     */
    xtype = ncadios_to_mpi_type(atype);
    if (xtype != itype){
        esize = adios_type_size (atype, adata);
        nelems = asize / esize;
        err = ncadiosiconvert(adata, buf, xtype, itype, nelems);
        if (err != NC_NOERR){
            return err;
        }
    }
    else{
        memcpy(buf, adata, asize);
    }

    free(adata);

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
    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);
}
