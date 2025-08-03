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
// #include "baseline_ncx_lib.h"
#include "baseline_ncx_app.h"
#include <malloc.h>
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



/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
#ifdef MEM_TRACKING
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
#endif

/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr(struct hdr *hdr_data, int ncid, int rank){
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int *dimid = (int *)malloc(ndims * sizeof(int));
    int i,j,k,nerrs=0;
    int err;

    for (i=0; i<ndims; i++){
        def_start_time = MPI_Wtime();
        err = ncmpi_def_dim(ncid, hdr_data->dims.value[i]->name,  hdr_data->dims.value[i]->size, &dimid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        // }
    }
    //pnetcdf_check_crt_mem(MPI_COMM_WORLD, -2);
    //define variables
    int nvars = hdr_data->vars.ndefined;
    int *varid = (int *)malloc(nvars * sizeof(int));
    int v_ndims, v_namelen, xtype, n_att;
    int *v_dimids;
    int att_namelen, att_xtype, att_nelems;

    for (i=0; i<nvars; i++){

        v_namelen =  hdr_data->vars.value[i]->name_len;
        xtype = hdr_data->vars.value[i]->xtype;

        v_ndims = hdr_data->vars.value[i]->ndims;
        v_dimids = (int *)malloc(v_ndims * sizeof(int));
        for(j=0; j<v_ndims; j++) v_dimids[j] = dimid[hdr_data->vars.value[i]->dimids[j]];
        def_start_time = MPI_Wtime();
        err = ncmpi_def_var(ncid, hdr_data->vars.value[i]->name, xtype, v_ndims,  v_dimids, &varid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        n_att = hdr_data->vars.value[i]->attrs.ndefined;
        // printf("\nn_att: %d\n", n_att);
        for(k=0; k<n_att; k++){
            att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
            att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
            att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
            err = ncmpi_put_att(ncid, varid[i],  hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
            att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
        }
        free(v_dimids);


    }
    free(varid);
    free(dimid);
    return nerrs;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size, status, err, nerrs=0;
    double end_time, start_time, start_time1, end_time1, start_time2, start_time3, max_time, min_time;
    double io_time, enddef_time, close_time, end_to_end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc < 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <source_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    source_name = argv[1];
    output_name = argv[2];

    if (rank == 0) {
        printf("Input file: %s\n", source_name);
        printf("Output file: %s\n", output_name);
    }
    struct hdr all_hdr;
    read_metadata_from_file(source_name, &all_hdr);
    struct hdr local_hdr;
    distribute_metadata(rank, size, &all_hdr, &local_hdr);
    free_hdr_meta(&all_hdr);

    // struct hdr recv_hdr;
    // create_dummy_data(rank, &dummy);


    

    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(source_name) - 3;
    // strncpy(filename, source_name, position);
    // strcat(filename, "_new");
    // strcat(filename, source_name + position);
    // if (rank==0) printf("\n%s\n", output_name);
    // print_memory_info();
    //app_check_crt_mem(MPI_COMM_WORLD, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // MPI_Info_set(info, "nc_hash_size_dim", "16777216");
    // MPI_Info_set(info, "nc_hash_size_var", "8388608");
    // MPI_Info_set(info, "nc_hash_size_dim", "4096");
    // MPI_Info_set(info, "nc_hash_size_var", "4096");
    
    err = ncmpi_create(MPI_COMM_WORLD, output_name, cmode, info, &ncid); ERR
    MPI_Barrier(MPI_COMM_WORLD);

    // printf("rank %d, recv_displs: %d, recvcounts: %d \n",  rank, recv_displs[i], recvcounts[i]);

    start_time1 = MPI_Wtime();
    define_hdr(&local_hdr, ncid, rank);
    io_time = MPI_Wtime() - start_time1;
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 1);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 1);
#endif
    free_hdr_meta(&local_hdr);
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 2);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 2);
#endif
    start_time2 = MPI_Wtime();
    err = ncmpi_enddef(ncid); ERR
    enddef_time = MPI_Wtime() - start_time2;
    
    // Clean up
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 3);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 3);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    start_time3 = MPI_Wtime();
    // print_memory_info();
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 4);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 4);
#endif
    err = ncmpi_close(ncid); ERR
    end_time =  MPI_Wtime();
    close_time = end_time - start_time3;
    end_to_end_time = end_time - start_time;
    
    // //pnetcdf_check_crt_mem(MPI_COMM_WORLD, 2);

    double times[5] = {end_to_end_time, io_time, enddef_time, close_time, total_def_time};
    char *names[5] = {"end-end", "create", "enddef (consistency check)", "close", "def_dim/var"};
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
    free_allocation_struct();
#endif
    
    MPI_Finalize();
    return 0;
}