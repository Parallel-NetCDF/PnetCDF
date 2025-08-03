#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include <inttypes.h>
#include <sys/resource.h>
#include "baseline_ncx_app.h"




#define FAIL -1

#define DIRTY_BYTES_THRESHOLD (256 * 1024) // 256kb

const char *src_file = NULL;
const char *out_file = NULL;
double crt_start_time,crt_time=0;
int group_create_count = 0;
int dataset_create_count = 0;

//used to avoid open & closing same group for every dataset

char last_group_name[256] = ""; 
// Define a structure to represent a dataset
// typedef struct {
//     int name_len;
//     int dtype_size;
//     int dspace_size;
//     char *name;
//     char *dtype;
//     char *dspace;
// } h5_dataset;

// // Define a structure to represent a group
// typedef struct {
//     char *name;
//     int name_len;
//     int ndst;
//     h5_dataset **datasets;
// } h5_group;

// // Define a structure to represent an array of groups
// typedef struct {
//     int ngrps;
//     h5_group **groups;
// } h5_grouparray;


// void print_datatype(hid_t datatype_id) {
//     if (H5Tget_class(datatype_id) == H5T_INTEGER) {
//         printf("%"PRId64 "is Integer\n",datatype_id);
//     } else if (H5Tget_class(datatype_id) == H5T_FLOAT) {
//         printf("%"PRId64 "is Float\n",datatype_id);
//     } else if (H5Tget_class(datatype_id) == H5T_STRING) {
//         printf("%"PRId64 "is String\n",datatype_id);
//     } else {
//         printf("%"PRId64 "is Other\n",datatype_id);
//     }
// }
// Function to create a group

// Helper to extract group and dataset names from full variable name
void split_var_name(const char* full_name, char* group_name, char* dataset_name) {
    const char* ptr = full_name;
    int underscore_count = 0;
    const char* split_ptr = NULL;

    while (*ptr != '\0') {
        if (*ptr == '_') {
            underscore_count++;
            if (underscore_count == 5) {
                split_ptr = ptr;
                break;
            }
        }
        ptr++;
    }

    if (split_ptr) {
        int group_len = split_ptr - full_name;
        strncpy(group_name, full_name, group_len);
        group_name[group_len] = '\0';

        // remove "_grp_", 

        memmove(group_name, group_name + 5, strlen(group_name + 5) + 1);

        strcpy(dataset_name, split_ptr + 1);
    } else {
        strcpy(group_name, "/");
        strcpy(dataset_name, full_name);
    }
}

static int deserialize_all_hdr(struct hdr **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (struct hdr *)malloc(sizeof(struct hdr));
        deserialize_hdr_meta(all_recv_hdr[i], all_collections_buffer + recv_displs[i], recvcounts[i]);
    }
    return 0;
}



int define_hdr_hdf5(struct hdr *hdr_data, hid_t file_id) {
    int nerrs = 0;
    int nvars = hdr_data->vars.ndefined;

    hid_t current_group_id = -1;
    // For tracking group context
    
    

    for (int i = 0; i < nvars; i++) {
        hdr_var *var = hdr_data->vars.value[i];
        char group_name[256], dataset_name[256];
        split_var_name(var->name, group_name, dataset_name);

        // Create or reuse group
        if (strcmp(group_name, last_group_name) != 0) {
            //need to create a new group, the group name has changed
            if (current_group_id >= 0)
                H5Gclose(current_group_id);
            crt_start_time = MPI_Wtime();
            current_group_id = H5Gcreate2(file_id, group_name,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            crt_time += MPI_Wtime() - crt_start_time;
            strcpy(last_group_name, group_name);
        } else if(i == 0){
            // If it's the first variable and the group is the same as the last one, we need to reopen the group
            current_group_id = H5Gopen2(file_id, group_name, H5P_DEFAULT);
        }
        // Map NC_ type to HDF5 native type
        hid_t h5type;
        switch (var->xtype) {
            case NC_INT:    h5type = H5T_NATIVE_INT; break;
            case NC_FLOAT:  h5type = H5T_NATIVE_FLOAT; break;
            case NC_DOUBLE: h5type = H5T_NATIVE_DOUBLE; break;
            case NC_CHAR:   h5type = H5T_NATIVE_CHAR; break;
            default:        h5type = H5T_NATIVE_INT; break;
        }

        // Create dataspace using dimension sizes
        hsize_t *dims = (hsize_t *)malloc(var->ndims * sizeof(hsize_t));
        for (int j = 0; j < var->ndims; j++) {
            dims[j] = hdr_data->dims.value[var->dimids[j]]->size;
        }
        hid_t space_id = H5Screate_simple(var->ndims, dims, NULL);
        free(dims);

        // Dataset creation property list
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_LATE);
        H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);


        


        // Create dataset in current group
        crt_start_time = MPI_Wtime();
        hid_t dset_id = H5Dcreate2(current_group_id, dataset_name, h5type, space_id,
                                   H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        crt_time += MPI_Wtime() - crt_start_time;
        H5Pclose(dcpl_id);
        H5Sclose(space_id);

        

        // Create attributes
        for (int k = 0; k < var->attrs.ndefined; k++) {
            hdr_attr *att = var->attrs.value[k];

            hid_t att_type;
            switch (att->xtype) {
                case NC_INT:    att_type = H5T_NATIVE_INT; break;
                case NC_FLOAT:  att_type = H5T_NATIVE_FLOAT; break;
                case NC_DOUBLE: att_type = H5T_NATIVE_DOUBLE; break;
                case NC_CHAR:   att_type = H5T_C_S1; break;
                default:        att_type = H5T_NATIVE_INT; break;
            }

            hsize_t att_dims = att->nelems;
            hid_t att_space = H5Screate_simple(1, &att_dims, NULL);
            hid_t attr_id = H5Acreate2(dset_id, att->name, att_type, att_space,
                                       H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr_id, att_type, att->xvalue);
            H5Sclose(att_space);
            H5Aclose(attr_id);
        }

        H5Dclose(dset_id);
    }
    H5Gclose(current_group_id);

    return nerrs;
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
    hid_t outfile_id, plist_id, fcpl_id;
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
    src_file = argv[1];
    out_file = argv[2];

    if (rank == 0) {
        printf("Input file: %s\n", src_file);
        printf("Output file: %s\n", out_file);
    }
    struct hdr all_hdr;
    read_metadata_from_file(src_file, &all_hdr);
    struct hdr local_hdr;
    distribute_metadata(rank, nproc, &all_hdr, &local_hdr);
    

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
    
    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(src_file) - 3;
    // strncpy(filename, src_file, position);
    // strcat(filename, "_new");
    // strcat(filename, src_file + position);
    // if (rank==0) printf("\n%s\n", out_file);
    struct hdr **all_recv_hdr = (struct hdr **)malloc(nproc * sizeof(struct hdr*));

    deserialize_all_hdr(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);
    
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    int block_size = 4 * 1024 * 1024;
    unsigned ik = 32;
    unsigned lk = 5;
    fcpl_id = H5Pcreate(H5P_FILE_CREATE);
    // H5Pset_istore_k(fcpl_id, 1024);
    H5Pset_sym_k(fcpl_id, ik, lk);
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    H5Pset_coll_metadata_write(plist_id, true);
    H5Pset_meta_block_size(plist_id, block_size);

    //Set dirty bytes threshold to 16 MB
    H5AC_cache_config_t mdc_config;
    // Initialize and retrieve current metadata cache configuration
    memset(&mdc_config, 0, sizeof(mdc_config));
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(plist_id, &mdc_config);

    // Set the dirty bytes threshold to 16 MB, so all processes will reach sync point until this amount of matadata cache is reached. Default: 256kb
    mdc_config.dirty_bytes_threshold = DIRTY_BYTES_THRESHOLD;
    // mdc_config.metadata_write_strategy = H5AC_METADATA_WRITE_STRATEGY__PROCESS_0_ONLY;

    // Apply the new configuration
    H5Pset_mdc_config(plist_id, &mdc_config);
    outfile_id = H5Fcreate(out_file, H5F_ACC_TRUNC, fcpl_id, plist_id);
    for (int i = 0; i < nproc; i++) {
        define_hdr_hdf5(all_recv_hdr[i], outfile_id);
    }
    
    io_time = MPI_Wtime() - end_time1;
    H5Pclose(plist_id);
    H5Pclose(fcpl_id);
    free_all_hdr(all_recv_hdr, nproc);
    MPI_Barrier(MPI_COMM_WORLD);
    
    end_time3 = MPI_Wtime();
    H5Fclose(outfile_id);
    end_time = MPI_Wtime();
    close_time = end_time - end_time3;
    end_to_end_time = end_time - start_time;
    mpi_time = end_time1 - start_time1;
    double times[5] = {end_to_end_time, mpi_time, io_time, close_time, crt_time};
    char *names[5] = {"end-end", "metadata exchange", "create (consistency check)", "close", "H5Dcreate & H5Gcreate"};
    double max_times[5], min_times[5];


    MPI_Reduce(&times[0], &max_times[0], 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) printf("ik: %u, lk: %u\n", ik, lk);
    for (int i = 0; i < 5; i++){
        if (rank == 0) {
            
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
    for (int i = 0; i < 5; i++){
        if (rank == 0) {
            printf("%f\n", names[i], max_times[i]);
            printf("%f\n", names[i], min_times[i]);
        }
    }
    //     if (rank == 0) {
    //     printf("H5Gcreate2 called %d times on rank 0\n", group_create_count);
    //     printf("H5Dcreate called %d times on rank 0\n", dataset_create_count);
    // }
    // Free memory used by the group array
    free_hdr_meta(&local_hdr);
    free_hdr_meta(&all_hdr);
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    free(recvcounts);
    MPI_Finalize();
    return 0;
}
