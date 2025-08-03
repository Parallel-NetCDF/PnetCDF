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
#include "baseline_ncx_app.h" 
#include <math.h>
#include <malloc.h>
#include "mem_tracker.h"

#ifdef MEM_TRACKING
#define malloc(size) tracked_malloc(size)
#define free(ptr)    tracked_free(ptr)
#define realloc(ptr, size) tracked_realloc(ptr, size)
#define strdup(s) tracked_strdup(s)
#endif

static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

const char *file_name = NULL;
const char *output_name = NULL;
double def_start_time, total_def_time=0;

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
int define_hdr(struct hdr *hdr_data, int ncid){
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

static int deserialize_all_hdr(struct hdr **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (struct hdr *)malloc(sizeof(struct hdr));
        deserialize_hdr_meta(all_recv_hdr[i], all_collections_buffer + recv_displs[i], recvcounts[i]);
    }
    return 0;
}

int define_all_hdr(struct hdr **all_recv_hdr, int nproc, int ncid){
    for (int i=0; i< nproc; i++){
        struct hdr *hdr_data = all_recv_hdr[i];
        define_hdr(hdr_data, ncid);
    }
    return 0;
}

static int free_all_hdr(struct hdr **all_recv_hdr, int nproc){
    if (all_recv_hdr != NULL){
        for (int i=0; i< nproc; i++){
            free_hdr_meta(all_recv_hdr[i]);
            free(all_recv_hdr[i]);
        }
        free(all_recv_hdr);
    }
    return 0;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, nerrs=0;
    double end_to_end_time, mpi_time, io_time, enddef_time, close_time, max_time, min_time;
    double start_time, start_time1, end_time1, end_time2, end_time3, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (argc < 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <source_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    file_name = argv[1];
    output_name = argv[2];

    if (rank == 0) {
        printf("Input file: %s\n", file_name);
        printf("Output file: %s\n", output_name);
    }
    struct hdr all_hdr;
    read_metadata_from_file(file_name, &all_hdr);
    struct hdr local_hdr;
    distribute_metadata(rank, nproc, &all_hdr, &local_hdr);
    free_hdr_meta(&all_hdr);
    // struct hdr recv_hdr;
    // create_local_hdr_data(rank, &local_hdr);

   // Print the created data for each process
    // printf("\nRank %d:\n", rank);
    // printf("Total Header Size: %lld\n", local_hdr.xsz);
    // printf("Dimensions:\n");
    // for (int i = 0; i < local_hdr.dims.ndefined; i++) {
    //     printf("  Name: %s, Size: %lld\n", local_hdr.dims.value[i]->name, local_hdr.dims.value[i]->size);
    // }

    // printf("Variables:\n");
    // for (int i = 0; i < local_hdr.vars.ndefined; i++) {
    //     printf("  Name: %s, Type: %d, NumDims: %d\n", local_hdr.vars.value[i]->name,  local_hdr.vars.value[i]->xtype, 
    //     local_hdr.vars.value[i]->ndims);
    //     printf("    Dim IDs: ");
    //     for (int j = 0; j < local_hdr.vars.value[i]->ndims; j++) {
    //         printf("%d ", local_hdr.vars.value[i]->dimids[j]);
    //     }
    //     printf("\n");
    //     printf("    Attributes:\n");
    //     for (int k = 0; k < local_hdr.vars.value[i]->attrs.ndefined; k++) {
    //         printf("      Name: %s, Nelems: %lld, Type: %d\n", local_hdr.vars.value[i]->attrs.value[k]->name, 
    //         local_hdr.vars.value[i]->attrs.value[k]->nelems, local_hdr.vars.value[i]->attrs.value[k]->xtype);
    //     }
    // }
    // printf("rank %d, buffer size: %lld \n", rank, local_hdr.xsz);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = start_time1 = MPI_Wtime();
    char* send_buffer = (char*) malloc(local_hdr.xsz);
    status = serialize_hdr_meta(&local_hdr, send_buffer);


    


    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) malloc(nproc * sizeof(MPI_Offset));
    MPI_Allgather(&local_hdr.xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

    // Calculate displacements for the second phase
    int* recv_displs = (int*) malloc(nproc * sizeof(int));
    int total_recv_size, min_size, max_size;
    total_recv_size = min_size = max_size = all_collection_sizes[0];
    recv_displs[0] = 0;


    for (int i = 1; i < nproc; ++i) {
        recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
        total_recv_size += all_collection_sizes[i];
        if(all_collection_sizes[i] > max_size){
            max_size = all_collection_sizes[i];
        }
        if(all_collection_sizes[i] < min_size){
            min_size = all_collection_sizes[i];
        }
    }
    double total_recv_size_MB = total_recv_size / (1024.0 * 1024.0);
    double min_size_MB = min_size / (1024.0 * 1024.0);
    double max_size_MB = max_size / (1024.0 * 1024.0);
    if(rank==0){
        printf("\nTotal buffer size: %f MB", total_recv_size_MB);
        printf("\nMax buffer size: %f MB", max_size_MB);
        printf("\nMin buffer size: %f MB \n", min_size_MB);
    }
    
    // printf("\nrank %d, local_hdr xsz %lld", rank, local_hdr.xsz);
    // Allocate buffer for receiving all header data
    char* all_collections_buffer = (char*) malloc(total_recv_size);
    int* recvcounts =  (int*)malloc(nproc * sizeof(int));
    for (int i = 0; i < nproc; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    // Phase 2: Communicate the actual header data
    // Before MPI_Allgatherv
    MPI_Allgatherv(send_buffer, local_hdr.xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);
    // Deserialize the received data and print if rank is 0
    free(send_buffer);
    
    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(file_name) - 3;
    // strncpy(filename, file_name, position);
    // strcat(filename, "_new");
    // strcat(filename, file_name + position);
    // if (rank==0) printf("\n%s\n", output_name);
    struct hdr **all_recv_hdr = (struct hdr **)malloc(nproc * sizeof(struct hdr*));

    deserialize_all_hdr(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // MPI_Info_set(info, "nc_hash_size_dim", "1048576");
    // MPI_Info_set(info, "nc_hash_size_var", "1048576");
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 0);
#endif
    err = ncmpi_create(MPI_COMM_WORLD, output_name, cmode, info, &ncid); ERR
    MPI_Info_free(&info);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    define_all_hdr(all_recv_hdr, nproc, ncid);
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 1);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 1);
#endif

    // for (int i = 0; i < nproc; ++i) {
    //     struct hdr *recv_hdr = (struct hdr *)malloc(sizeof(struct hdr)); 
    //     deserialize_hdr_meta(recv_hdr, all_collections_buffer + recv_displs[i], recvcounts[i]);
    //     define_hdr(recv_hdr, ncid, rank);
    //     free_hdr_meta(recv_hdr);
    // }
    // pnetcdf_check_crt_mem(MPI_COMM_WORLD, 0);
    io_time = MPI_Wtime() - end_time1;

#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 2);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 2);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    end_time2 = MPI_Wtime();

    err = ncmpi_enddef(ncid); ERR

    end_time3 = MPI_Wtime();
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 3);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 3);
#endif

    enddef_time = end_time3 - end_time2;


    // Clean up
    free(recvcounts);
    free_hdr_meta(&local_hdr);
    free_all_hdr(all_recv_hdr, nproc);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time3 = MPI_Wtime();
    // print_memory_info();
#ifdef MEM_TRACKING
    app_check_crt_mem(MPI_COMM_WORLD, 4);
    pnetcdf_check_crt_mem(MPI_COMM_WORLD, 4);
#endif
    err = ncmpi_close(ncid); ERR
    end_time = MPI_Wtime();
    close_time = end_time - end_time3;
    end_to_end_time = end_time - start_time;
    mpi_time = end_time1 - start_time1;


    double times[6] = {end_to_end_time, mpi_time, io_time, enddef_time, total_def_time, close_time};
    char *names[6] = {"end-end", "metadata exchange", "create (consistency check)", "enddef", "def_dim/var", "close"};
    double max_times[6], min_times[6];


    MPI_Reduce(&times[0], &max_times[0], 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 6, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 6; i++) {
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
    free_allocation_struct();
    MPI_Finalize();
    return 0;
}