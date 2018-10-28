/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include "netcdf.h"
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "ncadios_driver.h"
#define ERR(e){if(e){printf("Error:%d\n",e);return 2;}}
#define DIVIDER "\t---------------------------------\n"
//#define DIVIDER "\t************************************\n"

struct var_dim
{
    uint16_t id;
    uint64_t rank;
    int      nc_dimid;
    char dimname[256];
};

void copy_buffer(struct adios_bp_buffer_struct_v1 *dest
                ,struct adios_bp_buffer_struct_v1 *src) {

    memcpy (dest, src, sizeof(struct adios_bp_buffer_struct_v1));
}

static int verbose=0;
static int is_input_fortran = 0; // 0 = C generated BP file, 1 = Fortran generated BP file

int ncd_gen_name (char *fullname, char *path, char *name) {
    int i;
    char *new_path = strdup (path);
    if ( path[0] == '/')
         new_path=new_path+1;
          
    for ( i = 0; i < strlen (new_path); i++) {
        if ( new_path[i] == '[' || new_path[i] == ']' || new_path[i] == '/' || new_path[i] == '\\')
            new_path[i] = '_';
    }
    if (*new_path != '\0') {
        if (new_path[i-1]!='_') {
            if (strcmp(name,"") )  
                sprintf (fullname, "%s_%s", new_path, name);
            else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0'; 
            }
        }
        else 
            if (strcmp(name,"") )  
                sprintf (fullname, "%s%s", new_path, name);
            else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0'; 
            }
    }
    else
        strcpy (fullname, name);
    return 0;
}
int ncd_attr_str_ds (NC_ad* ncid
                    ,struct adios_attribute_struct_v1 * attribute
                    ,struct adios_bp_buffer_struct_v1 * ptr_buffer
                    ,int count
                    ,struct var_dim * var_dims
                    ,int var_dims_count)
{
    char fullname[255];
    char *path = attribute->path;
    char *name = attribute->name;
    int  valid,retval,attid;

    ncd_gen_name (fullname, path, name);
    valid = -1;
    if (strcmp(path,"/")==0) {
        valid = NC_GLOBAL;
        strcpy(fullname, name);
    }
    else {
        ncd_gen_name (fullname, path, "");
        //printf("\tpathname:%s\n",fullname);
        retval=ncadiosi_inq_varid(ncid,fullname,&valid);
        if(retval < 0)
           return 1; 
        else
            strcpy(fullname, name);
        //ERR(retval); 
        //if (valid<0)
        //    ncd_gen_name (fullname, path,name);
        //printf("\t attr_name:%s %d\n",fullname,valid);
    }
    //retval=ncadiosi_inq_attid(ncid,valid,fullname,&attid);
    //printf("\tretval:%d attid=%d\n",retval,attid);
    //printf(DIVIDER);
    if (retval == NC_NOERR ) {
       //printf("\tattribute (%s) existed\n", fullname);
       return 0;
     }
    //else
       //printf("\tattribute: %s \n", fullname);

    //ncadiosi_redef(ncid);

    void *value = attribute->value;
    size_t len = 1; 
    enum ADIOS_DATATYPES type =  attribute->type;
    struct adios_var_payload_struct_v1 var_payload;  
    struct adios_index_var_struct_v1 * vars_root = 0;
//    struct adios_index_attribute_struct_v1 * atts_root = 0;

    var_payload.payload = 0;

    if ( attribute->is_var == adios_flag_yes) {
        adios_posix_read_vars_index (ptr_buffer);
        adios_parse_vars_index_v1 (ptr_buffer, &vars_root, NULL, NULL);
        //print_vars_index ( vars_root);
        while (vars_root) {
            if (vars_root->id == attribute->var_id) {
                type = vars_root->type;
                if (!(vars_root->characteristics->dims.dims)) { 
                    value = vars_root->characteristics->value; 
                    //printf("\t      var: %s = ", vars_root->var_name);
                }
                else {
                	// this cannot happen because ADIOS does not let writing attributes pointing to non-scalar variable
                    //ncadiosi_enddef(ncid);
                    return 1;
                }
                break; 
            } 
            vars_root = vars_root->next;
        }
    }
        //printf("\t      XML: ");   
    switch (type) {
         case adios_unsigned_byte:
            //retval=ncadiosi_put_att_uchar(ncid,valid,fullname,NC_BYTE,len,value);
            break;
         case adios_byte:
            //retval=ncadiosi_put_att_schar(ncid,valid,fullname,NC_BYTE,len,value);
            break;
         case adios_string:
            //printf("%s\n", (char *) value);    
            //retval=ncadiosi_put_att_text(ncid,valid,fullname, strlen(value),value);
            break;
         case adios_short:
            //printf("\tvaule: %s\n", *(short *) value);    
            //retval=ncadiosi_put_att_short(ncid,valid,fullname,NC_SHORT,len,value);
            ERR(retval); 
            break;
         case adios_integer:
            //printf("%d\n", *((int *) value));    
            //retval=ncadiosi_put_att_int(ncid,valid,fullname,NC_INT,len,value);
            break;
         case adios_long:
            //printf("\tvaule: %s\n", *(long *) value);    
            //retval=ncadiosi_put_att_long(ncid,valid,fullname,NC_LONG,len,value);
            break;
         case adios_real:
            //printf("\tvaule: %s\n", *(float *) value);    
            //retval=ncadiosi_put_att_float(ncid,valid,fullname,NC_FLOAT,len,value);
            break;
         case adios_double:
            //printf("\tvaule: %s\n", *(double *) value);    
            //retval=ncadiosi_put_att_double(ncid,valid,fullname,NC_DOUBLE,len,value);
            break;
         default:
            break;
     }
     ERR(retval);

    //ncadiosi_enddef(ncid); 
    if ( var_payload.payload)
        free (var_payload.payload);
    return 0;
}

int ncd_dataset (NC_ad* ncid
                ,struct adios_var_header_struct_v1 *ptr_var_header
                ,struct adios_var_payload_struct_v1 *ptr_var_payload
                ,struct adios_bp_buffer_struct_v1 * ptr_buffer
                ,struct var_dim *var_dims
                ,int var_dims_count) {

    char *name = ptr_var_header->name;
    char *path = ptr_var_header->path;
    char fullname[256],dimname[256];
    enum ADIOS_DATATYPES type = ptr_var_header->type;
    //enum ADIOS_FLAG is_dim = ptr_var_header->is_dim;
    void *val = ptr_var_payload->payload; 
    //uint64_t payload_size = ptr_var_header->payload_size;
    struct adios_dimension_struct_v1 *dims = ptr_var_header->dims; 
    int maxrank = 0, i,j, valid=-1, nc_dimid=-1, retval=0;
    size_t rank = 0, start_dims[10],count_dims[10];
    int dimids[10];
    static int onename_dimid = -1;
    //int time_index, time_dimrank=0;
    struct adios_index_attribute_struct_v1 * atts_root = 0;
    
    memset(dimids,-1,10*sizeof(int));
    if(!strcmp(path,"") && !strcmp(name,""))
        return 0;
    ncd_gen_name (fullname, path, name);
   
    //printf(DIVIDER);
    //printf("\t  dataset: %s\n", fullname);

    val = ptr_var_payload->payload;
    //ncadiosi_redef(ncid);
    enum ADIOS_FLAG time_flag;
    while (dims) {
        ++maxrank;
        if (dims->dimension.is_time_index == adios_flag_yes) {
            time_flag = adios_flag_yes;
            //time_dimrank = maxrank-1;
            //fprintf(stderr, "Time dim rank = %d\n",time_dimrank); 
        } 
        dims = dims->next;
    } 
    //fprintf(stderr, "variable %s, ndims = %d time_flag = %d\n",ptr_var_header->name, maxrank, time_flag); 
    dims = ptr_var_header->dims;
    //time_index = 0;
    if (dims) {
        for (rank = 0; rank < maxrank; rank++) {
            //            printf("rank cal:time_index=%d j=%d rank=%d\n",time_dimrank,j,rank);
            /**********************************************************************
             * Process dataset which has global bounds with dynamic dimension value
             **********************************************************************/
            if ( dims->global_dimension.var_id != 0 ) {
                for (i = 0; i < var_dims_count; i++) {
                    if (var_dims [i].id == dims->global_dimension.var_id) {
                        dimids [rank] = var_dims [i]. nc_dimid;
                        break;
                    }
                }
                if (i==var_dims_count) {
                    adios_posix_read_attributes_index (ptr_buffer);
                    adios_parse_attributes_index_v1 (ptr_buffer, &atts_root);
                    while (atts_root) {
                        if (atts_root->id == dims->global_dimension.var_id) {
                            dimids[ rank] = *(int*)atts_root->characteristics->value;
                            break; 
                        }
                        atts_root = atts_root->next; 
                    }
                }
                if (dims->dimension.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->dimension.var_id) {
                            count_dims [ rank]=var_dims [i].rank;
                            break;
                        }
                    }
                    if (i==var_dims_count) {
                        adios_posix_read_attributes_index (ptr_buffer);
                        adios_parse_attributes_index_v1 (ptr_buffer, &atts_root);
                        //print_attributes_index ( atts_root);
                        while (atts_root) {
                            if (atts_root->id == dims->dimension.var_id) {
                                count_dims [ rank] = *(int*)atts_root->characteristics->value;
                                break; 
                            } 
                            atts_root = atts_root->next;
                        }
                    }
                }
                else
                    count_dims [ rank] = dims->dimension.rank;

                if ( dims->local_offset.var_id != 0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->local_offset.var_id){
                            start_dims[rank]=var_dims [i]. rank;;
                            break;
                        }
                    }
                    if (i==var_dims_count) {
                        adios_posix_read_attributes_index (ptr_buffer);
                        adios_parse_attributes_index_v1 (ptr_buffer, &atts_root);
                        while (atts_root) {
                            if (atts_root->id == dims->local_offset.var_id) {
                                start_dims [rank] = *(int*)atts_root->characteristics->value;
                                break; 
                            } 
                            atts_root = atts_root->next;
                        }
                    }
                }
                else{
                    start_dims[ rank]=dims->local_offset.rank;
                }
            }
            /**********************************************************************
             * Process dataset which has global bounds with constant dimension value
             ***********************************************************************/
            else if (dims->global_dimension.rank !=0 ) {
                //printf(" \tconstant global_info: %s rank: %d\n",fullname, dimids[rank]);
                dimids[rank] = dims->global_dimension.rank;
                if (dims->dimension.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->dimension.var_id)
                            count_dims[rank]=var_dims [i]. rank;
                    }
                }
                else
                    dimids[rank] = dims->dimension.rank;
                if (dims->local_offset.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->local_offset.var_id)
                            start_dims[rank]=var_dims [i]. rank;
                    }
                }
                else
                    start_dims[rank]=dims->local_offset.rank;
            }
            /*******************************************
             * Process dataset which has no global bounds
             ********************************************/
            else {
                //printf("\tdim: %s %d\n",atts_root->attr_name,dims->dimension.rank);
                if ( dims->dimension.var_id!=0
                        ||time_flag == adios_flag_yes) {
                    if (dims->dimension.rank!=0) {
                        sprintf(dimname,"%s_%zu",fullname,rank);
                        dimids[rank]=-1;
                        ncadiosi_inq_dimid(ncid, dimname, &dimids[rank]); 
                        if (dimids [rank] <= 0) 
                            retval=ncadiosi_def_dim (ncid, dimname,dims->dimension.rank,&dimids[rank]);
                        start_dims[rank] = 0;
                        count_dims[rank] = dims->dimension.rank;
                        /*fprintf(stderr,"\tdim[%zu]: c(%zu):s(%zu): dimid=%d\n"
                                ,rank 
                                ,count_dims[rank] 
                                ,start_dims[rank]
                                ,dimids[rank]
                               );*/

                    }
                    else {
                        for (i = 0; i < var_dims_count; i++) {
                            if (var_dims [i].id == dims->dimension.var_id) {
                                if (dims->dimension.is_time_index == adios_flag_yes) {
                                    start_dims[rank] = var_dims[i].rank - 1;
                                    //time_index = var_dims[i].rank;
                                    count_dims[rank] = 1;
                                    dimids[rank] = var_dims [i].nc_dimid; 
                                    /*printf("\tdim[%d]: c(%d):s(%d): dimid=%d (time-index)\n"
                                      ,rank
                                      ,count_dims[rank]
                                      ,start_dims[rank]
                                      ,dimids[rank]
                                      ); 
                                     */
                                }
                                else {
                                    start_dims[rank] = 0;
                                    count_dims[rank] = var_dims[i].rank;
                                    dimids[rank]=var_dims[i].nc_dimid;
                                    /*fprintf(stderr,"\tdim[%zu]: c(%zu):s(%zu): dimid=%d\n"
                                            ,rank
                                            ,count_dims[rank]
                                            ,start_dims[rank]
                                            ,dimids[rank]
                                           ); */
                                } 
                                break;
                            }
                        }
                    }
                    if (i==var_dims_count) {
                        adios_posix_read_attributes_index (ptr_buffer);
                        adios_parse_attributes_index_v1 (ptr_buffer, &atts_root);
                        //print_attributes_index(atts_root);
                        while (atts_root) {
                            if (atts_root->id == dims->dimension.var_id) {
                                ncd_gen_name (dimname, atts_root->attr_path
                                        ,atts_root->attr_name); 
                                if (!atts_root->characteristics->value) {
                                    for (i = 0; i < var_dims_count; i++) {
                                        if (var_dims [i].id == atts_root->characteristics->var_id) {
                                            start_dims[rank]=0;
                                            count_dims [rank] = var_dims [i].rank;

                                            ncadiosi_inq_dimid(ncid, dimname, &dimids[rank]);
                                            if (dimids [rank] <= 0)
                                                ncadiosi_def_dim (ncid, dimname
                                                        ,var_dims[i].rank
                                                        ,&dimids [rank]);
                                            i = var_dims_count + 1; 
                                        }                
                                    }                
                                }                
                                else {
                                    count_dims [ rank] = *(int *)atts_root->characteristics->value;
                                    start_dims[rank]=0;
                                    ncadiosi_inq_dimid(ncid, dimname, &dimids[rank]);
                                    if (dimids [rank] <= 0)
                                        ncadiosi_def_dim (ncid, dimname
                                                ,*(int *)atts_root->characteristics->value
                                                ,&dimids [rank]); 
                                }
                                /*printf("\t local[%d]: c(%d) id(%d)\n"
                                  ,rank,count_dims[rank], dimids[rank]); */
                                break; 
                            } 
                            atts_root = atts_root->next;
                        }
                    }
                }
                else {
                    //printf ("Error, every dimension in netcdf need to have name!\n");
                    sprintf(dimname,"%s_%zu", fullname,rank);
                    ncadiosi_inq_dimid(ncid,dimname,&nc_dimid);
                    if (nc_dimid<0)
                        retval = ncadiosi_def_dim ( ncid, dimname, dims->dimension.rank, &nc_dimid);
                    dimids[rank]=nc_dimid;
                    count_dims[rank] = dims->dimension.rank;
                    start_dims[rank] =0; 
                    ERR(retval);
                    fprintf(stderr, "\t local[%zu]: %" PRIu64 "\n",rank,dims->dimension.rank);
                } 
            }
            if (dims)
                dims = dims->next;
        } // end of for loop
        val = ptr_var_payload->payload;    
        ncadiosi_inq_varid(ncid,fullname,&valid);
        //ncadiosi_redef(ncid);
        int time_idx=-1;
        for (rank = 0; rank < maxrank; rank++) {
            if (count_dims[rank]==0) {
                count_dims[rank]=1;
                time_idx = rank;
                break;
            }
        } 
        for (rank = 0; rank < maxrank; rank++) {
            if (verbose>0)
                fprintf(stderr, "\tdimension info[%zu]: c(%zu) s(%zu)\n"
                        ,rank,count_dims[rank], start_dims[rank]);
        }
        if (time_idx == 0 && dimids[time_idx]!=0) {
            for (rank=maxrank-1;rank>0;rank--) {
                dimids[rank] = dimids[rank-1];                
                start_dims[rank] = start_dims[rank-1];                
            }
            start_dims[time_idx] = var_dims[0].rank-1;
            dimids[time_idx] = var_dims[0].id;                
        }

        /* In case of Fortran generated file, we have to flip the dimensions here
           because the dimension order at this moment is how it was in Fortran, but
           this is a C code and thus NetCDF expects C ordering of dimensions here.
        */
        if (is_input_fortran==1) {
            size_t tmp;
            //fprintf(stderr, "\tFlip %d Fortran dimensions, time_idx=%d\n", maxrank, time_idx);
            for (rank=0; rank < maxrank; rank++) {
                if (rank < maxrank-1-rank) {
                    //fprintf(stderr, "\t  flip rank %d with rank %d\n", rank, maxrank-1-rank);
                    tmp = dimids[rank];
                    dimids[rank] = dimids[maxrank-1-rank];
                    dimids[maxrank-1-rank] = tmp;

                    tmp = count_dims[rank];
                    count_dims[rank] = count_dims[maxrank-1-rank];
                    count_dims[maxrank-1-rank] = tmp;

                    tmp = start_dims[rank];
                    start_dims[rank] = start_dims[maxrank-1-rank];
                    start_dims[maxrank-1-rank] = tmp;
                }
            }
        }

        switch(type) { 
            case adios_real:
                //for (rank=0;rank<maxrank;rank++)
                //         printf("\t local[%d]: c(%d) (%d) %d\n",
                //    rank,count_dims[rank], start_dims[rank], dimids[rank]);
                if ( valid<0) {
                    retval=ncadiosi_def_var(ncid,fullname,NC_FLOAT,maxrank,dimids,&valid);
                    ERR(retval);
                }
               // retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_vara_float(ncid,valid,start_dims,count_dims,val);
                ERR(retval);
                break;
            case adios_double:
                if ( valid<0) 
                    retval=ncadiosi_def_var(ncid,fullname,NC_DOUBLE,maxrank,dimids,&valid);
                ERR(retval);
                //retval=ncadiosi_enddef(ncid);
                //ERR(retval);
                //retval=ncadiosi_put_vara_double(ncid,valid,start_dims,count_dims,val);
                //ERR(retval);
                //printf("end   writing!\n");
                break;
            case adios_long:
                if ( valid<0) 
                    retval=ncadiosi_def_var(ncid,fullname,NC_LONG,maxrank,dimids,&valid);
                //retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_vara_long(ncid,valid,start_dims,count_dims,val);
                break;
            case adios_unsigned_byte:
                if ( valid<0) 
                    retval=ncadiosi_def_var(ncid,fullname,NC_BYTE,maxrank,dimids,&valid);
                //retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_vara_uchar(ncid,valid,start_dims,count_dims,val);
                break;
            case adios_byte:
                //printf("write byte test %d %d\n",maxrank,dimids[0]);
                if ( valid<0) 
                    retval=ncadiosi_def_var(ncid,fullname,NC_BYTE,maxrank,dimids,&valid);
                ERR (retval);
                //printf("\t vid=%d\n",valid);
                //retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_vara_schar(ncid,valid,start_dims,count_dims,val);
                //ERR (retval);
                //printf("write byte test\n");
                break;
            case adios_integer:
                if (valid < 0) {
                    retval = ncadiosi_def_var (ncid,fullname,NC_INT,maxrank,dimids,&valid);
                    //printf("definition done!\n");
                } 
                //retval = ncadiosi_enddef (ncid);
                //ERR (retval);
                //retval=ncadiosi_put_vara_int (ncid,valid,start_dims,count_dims,val);
                //ERR (retval);
                break;
            default:
                //retval=ncadiosi_enddef(ncid);
                break;
        }
    }
    else if (ptr_var_header->is_dim == adios_flag_yes) {
        for ( j = 0; j<var_dims_count;j++){
            if (var_dims [j].id==ptr_var_header->id) { 
                break; 
            }
        }

        //ncadiosi_redef(ncid);
        ncadiosi_inq_dimid ( ncid, fullname, &nc_dimid);
        if ( var_dims[j].rank == 0)
            return 0;
        if ( nc_dimid < 0) {
            retval = ncadiosi_def_dim ( ncid, fullname, var_dims[j].rank, &nc_dimid);
            ERR(retval);
        }
        var_dims [j].nc_dimid = nc_dimid;
    }
    else {
        rank = 1;
        //ncadiosi_redef(ncid);
        if (onename_dimid==-1)
        {
            retval=ncadiosi_def_dim (ncid, "one", 1, &onename_dimid);
            ERR(retval);
        }
        else {
            //ncadiosi_redef(ncid);
            ncadiosi_inq_varid (ncid, fullname, &valid);
        }
        dimids[0]=onename_dimid;
        rank = 1;
        switch (type) {
            case adios_real:
                if (valid < 0 ) {
                    //printf("\t ncd-scalar-real: %d %d %s\n",dimids[0],valid, fullname);
                    retval=ncadiosi_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
                    ERR(retval);
                }
                //retval=ncadiosi_enddef(ncid);
                //ERR(retval);
                //retval=ncadiosi_put_var_float(ncid,valid,val);
                break;
            case adios_double:
                if (valid < 0 ) {
                    //printf("\t ncd-scalar: %d %d %s\n",dimids[0],valid, fullname);
                    retval=ncadiosi_def_var(ncid,fullname,NC_DOUBLE,rank,dimids,&valid);
                    ERR(retval);
                    //retval=ncadiosi_enddef(ncid);
                    //ERR(retval);
                }
                //retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_var_double(ncid,valid,val);
                break;
            case adios_long:
                if (valid < 0 ) {
                    //printf("\t ncd-scalar: %d %d %s\n",dimids[0],valid, fullname);
                    retval=ncadiosi_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
                    ERR(retval);
                }
                //retval=ncadiosi_enddef(ncid);
                //ERR(retval);
                retval=ncadiosi_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
                //retval=ncadiosi_enddef(ncid);
                //retval=ncadiosi_put_var_long(ncid,valid,val);
                break;
            case adios_integer:
                if (valid < 0 ) {
                    retval=ncadiosi_def_var(ncid,fullname,NC_INT,rank,dimids,&valid);
                    ERR(retval);
                }
                //retval=ncadiosi_enddef(ncid);
                //ERR(retval);
                //printf("\t   scalar: %d\n", *(int *)val);
                //retval=ncadiosi_put_var_int(ncid,valid,val);
                //ERR(retval);
                break;
            default:
                //retval=ncadiosi_enddef(ncid);
                break;
        }
    }
    return 0;
}

int ncadiosi_parse_header (NC_ad *ncid)
{
    char out_fname [256];
    int rc = 0;
    /*
    if (argc < 2)
    {
        fprintf (stderr, "usage: %s <argv[1]_in> [argv[1]_out]\n"
                ,argv [0]
                );

        return -1;
    }

    if (argc > 2)
        strcpy (out_fname,argv[2]);
    else 
    {
        int size = strlen(argv[1]);
        strcpy (out_fname,argv[1]);
        if(!strcmp (&out_fname[size-3], ".bp")) {
          out_fname [size-2] = 'n'; 
          out_fname [size-1] = 'c';
        }
        else {
          strcat(out_fname, ".nc");
        }
        
    }
    */
    //ncadiosi_create ( out_fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);

    struct adios_bp_buffer_struct_v1 * b = 0;
    struct adios_bp_buffer_struct_v1 * b_0 = 0;
    struct adios_bp_buffer_struct_v1 * b_1 = 0;
    uint32_t version = 0;
    is_input_fortran = 0;

    b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    b_0 = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    b_1 = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (b);

    rc = adios_posix_open_read_internal (ncid->path, "", b);
    if (!rc)
    {
        //fprintf (stderr, "bp2ncd: file not found: %s\n", argv[1]);

        return -1;
    }

    adios_posix_read_version (b);
    adios_parse_version (b, &version);

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;
    struct adios_index_attribute_struct_v1 * attrs_root = 0;

    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root, NULL);

    copy_buffer(b_0, b);
    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root, NULL, NULL);

    copy_buffer(b_1, b);
    adios_posix_read_attributes_index (b);
    adios_parse_attributes_index_v1 (b, &attrs_root);

    pg = pg_root;
    while (pg)
    {
        int i,j;
        int var_dims_count = 0;
        struct var_dim * var_dims = 0;

        struct adios_process_group_header_struct_v1 pg_header;
        struct adios_vars_header_struct_v1 vars_header;
        struct adios_attributes_header_struct_v1 attrs_header;

        struct adios_var_header_struct_v1 var_header;
        struct adios_var_payload_struct_v1 var_payload;
        struct adios_attribute_struct_v1 attribute;

        // setup here to read the process group from (and size)
        b->read_pg_offset = pg->offset_in_file;
        if (pg->next)
        {
            b->read_pg_size =   pg->next->offset_in_file
                              - pg->offset_in_file;
        }
        else
        {
            b->read_pg_size =   b->pg_index_offset
                              - pg->offset_in_file;
        }

        adios_posix_read_process_group (b);
        adios_parse_process_group_header_v1 (b, &pg_header);
        //printf ("*************************************************\n"); 
        //printf ("\tTime Index Name: %s %d\n", pg_header.time_index_name, pg_header.time_index);
        //printf ("*************************************************\n"); 

        // Note here if we deal with Fortran generated file
        if (pg_header.host_language_fortran == adios_flag_yes) {
            is_input_fortran = 1;
        }

        /****************************************
        * Create unlimited time index dimension 
        ****************************************/
        if (pg_header.time_index_name) {
             
             var_dims = realloc (var_dims, (var_dims_count + 1)
                          * sizeof (struct var_dim)
                          );
             static int time_dimid = -1;
             ncadiosi_def_dim(ncid,pg_header.time_index_name,NC_UNLIMITED,&time_dimid);
             //ncadiosi_enddef(ncid);
             strcpy(var_dims[var_dims_count].dimname,pg_header.time_index_name);
             var_dims[var_dims_count].id = 0; 
             var_dims[var_dims_count].rank = pg_header.time_index;
             var_dims[var_dims_count].nc_dimid = time_dimid;

             var_dims_count=var_dims_count+1;
        }

        adios_parse_vars_header_v1 (b, &vars_header);
  
        //printf("time-index id: %s %d\n",pg_header.time_index_name, vars_header.count);
        for (i = 0; i < vars_header.count; i++) {
            var_payload.payload = 0;
            adios_parse_var_data_header_v1 (b, &var_header);

            if (var_header.is_dim == adios_flag_yes) {
                var_payload.payload = malloc (var_header.payload_size);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                                ,var_header.payload_size
                                                );
            }
            else {
                // alloc size +1 to remove valgrind complaint
                var_payload.payload = malloc (var_header.payload_size + 1);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                                ,var_header.payload_size
                                                );
                ncd_dataset(ncid,&var_header, &var_payload,b_1,var_dims,var_dims_count);
            }

            if (var_header.is_dim == adios_flag_yes) {
                int flag=0;
                var_dims = realloc (var_dims,   (var_dims_count + 1)
                                              * sizeof (struct var_dim)
                                   );
                
                strcpy(var_dims [var_dims_count].dimname,var_header.name);
                var_dims [var_dims_count].id = -1;
                var_dims [var_dims_count].rank = -1;
                var_dims [var_dims_count].nc_dimid = -1;

                for( j = 0 ; j < var_dims_count; j++){
                    if (var_dims [j].id == var_header.id) { 
                         var_dims [j].rank = *(unsigned int *) var_payload.payload;
                         flag = 1;
                         break; 
                    }
                }

                if(flag == 0) {
                    var_dims [var_dims_count].id = var_header.id;
                    var_dims [var_dims_count].rank = *(unsigned int *)
                                                            var_payload.payload;
                    var_dims_count++;
                }
                ncd_dataset(ncid,&var_header, &var_payload,b_1,var_dims,var_dims_count);
            }

            if (var_payload.payload)
            {
                free (var_payload.payload);
            }
            //printf ("\n");
        }

        adios_parse_attributes_header_v1 (b, &attrs_header);

        for (i = 0; i < attrs_header.count; i++)
        {

            //adios_parse_attribute_v1 (b, &attribute);
            //ncd_attr_str_ds (ncid, &attribute, b_0, vars_header.count, var_dims, var_dims_count);
        }

        var_dims_count = 0;
        if (var_dims)
            free (var_dims);
        pg = pg->next;
    }
    //printf (DIVIDER);
    //printf ("End of %s\n", argv[1]);

    adios_posix_close_internal (b);
    free (b);
    free (b_0);
    //ncadiosi_close (ncid);
    return 0;
}

