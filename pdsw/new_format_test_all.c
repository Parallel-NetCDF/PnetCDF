/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>
#include <assert.h>
// #include "baseline_ncx_lib.h"
#include "baseline_ncx_app.h" 
#include "mem_tracker.h"

#ifdef MEM_TRACKING
#define malloc(size) tracked_malloc(size)
#define free(ptr)    tracked_free(ptr)
#define realloc(ptr, size) tracked_realloc(ptr, size)
#endif


static int verbose;
const char *source_name = NULL;
const char *output_name = NULL;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}


double def_start_time;
double total_def_time = 0;



#ifdef MEM_TRACKING
/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("total maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   sum_size, (float)sum_size /1048576);
            // printf("rank 0 maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
            //        malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            // printf("rank 1 maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
            //        malloc_size, (float)malloc_size /1048576);
        }
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_crt_mem(MPI_Comm comm, int checkpoint)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("checkpoint %d: total current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   checkpoint, sum_size, (float)sum_size /1048576);
            // printf("checkpoint %d: rank 0 current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
            //        checkpoint, malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            // printf("checkpoint %d: rank 1 current heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
            //        checkpoint, malloc_size, (float)malloc_size /1048576);
        }
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

/* check PnetCDF library internal memory usage */
static int
app_check_crt_mem(MPI_Comm comm, int checkpoint)
{
    int err=0, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    malloc_size = inq_malloc_use();
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("checkpoint %d: total current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   checkpoint, sum_size, (float)sum_size /1048576);
            // printf("checkpoint %d: rank 0 current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
            //        checkpoint, malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            // printf("checkpoint %d: rank 1 current heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
            //        checkpoint, malloc_size, (float)malloc_size /1048576);
        }
    }
    return nerrs;
}
static int
app_check_mem_usage(MPI_Comm comm)
{
    int err=0, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    malloc_size = inq_max_malloc_use();
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0){
            printf("total maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
                   sum_size, (float)sum_size /1048576);
            // printf("rank 0 maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
            //        malloc_size, (float)malloc_size /1048576);
        }else if (rank == 1){
            // printf("rank 1 maximum heap memory allocated by App internally is %lld bytes (%.2f MB)\n",
            //        malloc_size, (float)malloc_size /1048576);
        }
    }
    return nerrs;
}
#endif

/* ---------------------------------- Decode Metadata ----------------------------------------*/
// int define_hdr_nf_distribute(struct hdr *hdr_data, int ncid, int rank, int var_start, int var_count) {
//     //define dimensions
//     int ndims= hdr_data->dims.ndefined;
    
//     int i,j,k,nerrs=0;
//     int err;
//     def_start_time = MPI_Wtime();
//     char str_blk[20];
//     sprintf(str_blk, "blk_rank_%d", rank);
//     int blkid;

//     def_start_time = MPI_Wtime();
//     err = ncmpi_def_block(ncid, str_blk, &blkid);
//     total_def_time += MPI_Wtime() - def_start_time;

//     //define variables
//     int *varid = (int *)malloc(var_count * sizeof(int));
//     int v_ndims, v_namelen, xtype, n_att;
//     int att_namelen, att_xtype, att_nelems;

//     for (i=var_start; i<var_start+var_count; i++){
        
//         int varindex = i - var_start;
//         v_namelen =  hdr_data->vars.value[i]->name_len;
//         xtype = hdr_data->vars.value[i]->xtype;

//         v_ndims = hdr_data->vars.value[i]->ndims;
//         int *v_dimids = (int *)malloc(v_ndims * sizeof(int));
//         for (j=0; j<v_ndims; j++){
//             int dimindex = hdr_data->vars.value[i]->dimids[j];
//             def_start_time = MPI_Wtime();
//             err = ncmpi_def_dim(ncid, blkid, hdr_data->dims.value[dimindex]->name,  hdr_data->dims.value[dimindex]->size, &v_dimids[j]); ERR
//             total_def_time += MPI_Wtime() - def_start_time;
//             // }
//         }

//         def_start_time = MPI_Wtime();

//         if (v_ndims == 0)
//             err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, NULL, &varid[varindex]);  
//         else
//             err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, v_dimids, &varid[varindex]); 
//         ERR

        
//         total_def_time += MPI_Wtime() - def_start_time;
//         n_att = hdr_data->vars.value[i]->attrs.ndefined;
//         assert(n_att == 0);
//         // printf("\nn_att: %d\n", n_att);
//         // for(k=0; k<n_att; k++){
//         //     att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
//         //     att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
//         //     att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
//         //     err = ncmpi_put_att(ncid, varid[i], hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
//         //     att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
//         // }
//         free(v_dimids);

//     }
//     free(varid);
//     return nerrs;
// }


int define_hdr_nf(struct hdr *hdr_data, int ncid, int rank) {
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int var_start = 0;
    int var_count = hdr_data->vars.ndefined;
    int i,j,k,nerrs=0;
    int err;
    def_start_time = MPI_Wtime();
    char str_blk[20];
    sprintf(str_blk, "blk_rank_%d", rank);
    int blkid;

    def_start_time = MPI_Wtime();
    err = ncmpi_def_block(ncid, str_blk, &blkid);
    total_def_time += MPI_Wtime() - def_start_time;

    //define variables
    int *varid = (int *)malloc(var_count * sizeof(int));
    int v_ndims, v_namelen, xtype, n_att;
    int att_namelen, att_xtype, att_nelems;

    for (i=var_start; i<var_start+var_count; i++){
        
        int varindex = i - var_start;
        v_namelen =  hdr_data->vars.value[i]->name_len;
        xtype = hdr_data->vars.value[i]->xtype;

        v_ndims = hdr_data->vars.value[i]->ndims;
        int *v_dimids = (int *)malloc(v_ndims * sizeof(int));
        for (j=0; j<v_ndims; j++){
            int dimindex = hdr_data->vars.value[i]->dimids[j];
            def_start_time = MPI_Wtime();
            err = ncmpi_def_dim(ncid, blkid, hdr_data->dims.value[dimindex]->name,  hdr_data->dims.value[dimindex]->size, &v_dimids[j]); ERR
            total_def_time += MPI_Wtime() - def_start_time;
            // }
        }

        def_start_time = MPI_Wtime();

        if (v_ndims == 0)
            err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, NULL, &varid[varindex]);  
        else
            err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, v_dimids, &varid[varindex]); 
        ERR

        
        total_def_time += MPI_Wtime() - def_start_time;
        n_att = hdr_data->vars.value[i]->attrs.ndefined;
        assert(n_att == 0);
        // printf("\nn_att: %d\n", n_att);
        // for(k=0; k<n_att; k++){
        //     att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
        //     att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
        //     att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
        //     err = ncmpi_put_att(ncid, varid[i], hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
        //     att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
        // }
        free(v_dimids);

    }
    free(varid);
    return nerrs;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, nerrs=0;
    double end_time, start_time, start_time1, end_time1, start_time2, start_time3, max_time, min_time;
    double enddef_time, close_time, end_to_end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (argc < 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <source_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    source_name = argv[1];
    output_name = argv[2];
    struct hdr all_hdr;
    // struct hdr recv_hdr;
    // create_dummy_data(rank, &dummy);

    //the following function's memory usage should be selectively tracked - production application dont need source file for metadata
    read_metadata_from_file(source_name, &all_hdr);
    struct hdr local_hdr;
    distribute_metadata(rank, nproc, &all_hdr, &local_hdr);
    free_hdr_meta(&all_hdr);

    

    int ncid, cmode;
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(source_name) - 3;
    // strncpy(filename, source_name, position);
    // strcat(filename, "_new");
    // strcat(filename, source_name + position);
    // if (rank==0) printf("\n%s\n", output_name);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // MPI_Info_set(info, "nc_hash_size_dim", "16777216");
    // MPI_Info_set(info, "nc_hash_size_var", "8388608");
    // MPI_Info_set(info, "nc_hash_size_dim", "16384");
    // MPI_Info_set(info, "nc_hash_size_var", "16384");
    // The following memory footprint should be selectively tracked because only a fraction of metadata is used in actual application
    // For simplicity, we only record 1/nproc fraction of the reported mem usage 
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 0);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 0);
#endif
    err = ncmpi_create(MPI_COMM_WORLD, output_name, cmode, info, &ncid); ERR
    double create_time = MPI_Wtime() - start_time;
    // MPI_Barrier(MPI_COMM_WORLD);

    // printf("rank %d, recv_displs: %d, recvcounts: %d \n",  rank, recv_displs[i], recvcounts[i]);
    int nvars = local_hdr.vars.ndefined;
    // int nvars = 10;
    // int vars_per_process = nvars / nproc;
    // int remainder = nvars % nproc;
    // int start = rank * vars_per_process + (rank < remainder ? rank : remainder);
    // int count = vars_per_process + (rank < remainder ? 1 : 0);
    // if (rank == 0) 
    //     printf("\ntotal var: %d", nvars);
    // printf("\nrank %d, start %d, count %d\n", rank, start, count);
    start_time1 = MPI_Wtime();
    define_hdr_nf(&local_hdr, ncid, rank);
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 1);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 1);
#endif
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 2);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 2);
#endif
    // MPI_Barrier(MPI_COMM_WORLD);
    start_time2 = MPI_Wtime();
    err = ncmpi_enddef(ncid); ERR
    enddef_time = MPI_Wtime() - start_time2;
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 3);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 3);
#endif

    // Clean up
    // MPI_Barrier(MPI_COMM_WORLD);
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 4);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 4);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    start_time3 = MPI_Wtime();
    err = ncmpi_close(ncid); ERR
    end_time =  MPI_Wtime();
    close_time = end_time - start_time3;
    end_to_end_time = end_time - start_time;

    free_hdr_meta(&local_hdr);

    double times[5] = {end_to_end_time, create_time,  enddef_time, close_time, total_def_time};
    char *names[5] = {"end-end", "create (consistency check)", "enddef", "close", "def_dim/var"};
    double max_times[5], min_times[5];

    MPI_Reduce(&times[0], &max_times[0], 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 5; i++) {
        if (rank == 0) {
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 5);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 5);
    pnetcdf_check_mem_usage(MPI_COMM_WORLD);
    app_check_mem_usage(MPI_COMM_WORLD);
#endif
    
    MPI_Finalize();
    return 0;



}