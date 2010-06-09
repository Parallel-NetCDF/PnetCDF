/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <string.h>
#include <unistd.h>

/* need this for SIZEOF_INT so we get the correct printf output */
#include "ncconfig.h" 

#ifdef SIZEOF_INT
# if SIZEOF_INT == 4
#  define lld(x) (x)
# elif  SIZEOF_INT == 8
#  define lld(x) (long long)(x)
# endif
#endif

/* The file name is taken as a command-line argument. */

/* Measures the I/O bandwidth for writing/reading a 3D
   block-distributed array to a file corresponding to the global array
   in row-major (C) order.
   Note that the file access pattern is noncontiguous.
  
   Array size 128^3. For other array sizes, change array_of_gsizes below.*/
#define TEST_HANDLE_ERR(status)                                 \
{                                                               \
  if ((status) != NC_NOERR)                                     \
    printf( "%s\n", ncmpi_strerror((status)) );                 \
}


int main(int argc, char **argv)
{
    int i, j, m, array_of_gsizes[3],array_of_distribs[3];
    int order, nprocs, len, **buf, mynod;
    MPI_Offset bufcount;
    int array_of_dargs[3], array_of_psizes[3];
    int status;
    MPI_Offset sizes[3], array_of_starts[3], stride[3];
    char basename[50], filename[100];
    char dimname[20], varname[20];
    int ncid, dimids0[3], dimids1[3], rank_dim[3], *varid;
    MPI_Info info;
    MPI_Offset **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    int ndims = 3;
    int nvars = 10;
    int k, k_loop;
    MPI_Datatype *datatype_list;
    int length;
    int mvar_flag = 0;
    int *array_of_requests;
    int unlimit_flag;
    int *array_of_statuses;
    int buf_var[32] ={1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2,
                      3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4};

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynod);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    length = 2;
    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    nvars = 4;
//    strcpy(basename,  "/pvfs2/kgao/test_non_blocking/test");
    strcpy(basename,  "./test");

    order = MPI_ORDER_C;
    
    buf = (int **)malloc(nvars*sizeof(int*));
    if (buf == NULL){
	printf("buf malloc error\n");
        return 0;		
    }	
    varid = (int *)malloc(nvars*sizeof(int));
    if (varid == NULL){
	printf("varid malloc error\n");
        return 0;		
    }	
    bufcount_list = (MPI_Offset *)malloc(nvars*sizeof(MPI_Offset));
    if (bufcount_list == NULL){
	printf("bufcount_list malloc error\n");
        return 0;		
    }	
    starts_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (starts_list== NULL){
	printf("starts_list malloc error\n");
        return 0;		
    }	
    count_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (count_list == NULL){
	printf("count_list malloc error\n");
        return 0;		
    }	
    datatype_list = (MPI_Datatype*)malloc(nvars*sizeof(MPI_Datatype));
    if (datatype_list == NULL){
	printf("count_list malloc error\n");
        return 0;		
    }	
    
    array_of_requests = (int *)malloc(nvars*sizeof(int));
    array_of_statuses = (int *)malloc(nvars*sizeof(int));
    
    for (i=0; i<nvars; i++) {
        starts_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
        if (starts_list[i] == NULL){
	    printf("starts_list[%d] malloc error\n", i);
            return 0;		
        }	
        count_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
        if (count_list[i] == NULL){
            printf("count_list[%d] malloc error\n", i);
            return 0;		
        }	
    }
  
    array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;

    array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
    
    bufcount = 1;
    for (i=0; i<ndims; i++) {
        array_of_psizes[i] = 0;
        sizes[i] = length;
        bufcount *= length;
    }
    MPI_Dims_create(nprocs, ndims, array_of_psizes);

    /* subarray in each process is len x len x len */
    for (i=0; i<ndims; i++)
        array_of_gsizes[i] = length * array_of_psizes[i];

    /* mynd's process rank in each dimension (in MPI_ORDER_C) */
    rank_dim[2] =  mynod %  array_of_psizes[2];
    rank_dim[1] = (mynod /  array_of_psizes[2]) % array_of_psizes[1];
    rank_dim[0] =  mynod / (array_of_psizes[2]  * array_of_psizes[1]);

    /* starting coordinates of the subarray in each dimension */
    for (i=0; i<ndims; i++)
        array_of_starts[i] = length * rank_dim[i];

    for (i=0; i<nvars; i++) {
        for (j=0; j<ndims; j++) {
           starts_list[i][j] = array_of_starts[j];
           count_list[i][j]  = length;
        }
        bufcount_list[i] = bufcount;
        datatype_list[i] = MPI_INT;
    }

    for (i=0; i<nvars;i++){
        buf[i] = (int *) malloc(bufcount * sizeof(int));
        if (buf[i] == NULL){
            printf("buf[i]malloc error\n");
            return 0;		
        }	
	for (j=0; j<bufcount; j++)
		buf[i][j]=mynod+1;
    }

    MPI_Info_create(&info);

/*
 *  MPI_Info_set(info, "group_cyclic_fd", "enable");
 *  MPI_Info_set(info, "cb_buffer_size", "1024");
 *  MPI_Info_set(info, "cb_buffer_size", "16777216");
 *  MPI_Info_set(info, "romio_no_indep_rw", "true");
 *   MPI_Info_set(info, "romio_cb_write", "true");
 */

    for (k=0; k<=9; k++){
//      sprintf(filename, "nc_%d.%d.%d.nc", length, nvars, k);
        sprintf(filename, "%s.%d.%d.%d.nc", basename, length, nvars, k);
        status = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_OFFSET,
                        info, &ncid);
        TEST_HANDLE_ERR(status);
      /* define dimensions */
        for (i=0; i<ndims; i++){
            sprintf(dimname, "dim0_%d", i);
            ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids0[i]);
        }
        sprintf(dimname, "dim1_%d", 0);
        ncmpi_def_dim(ncid, dimname, NC_UNLIMITED, &dimids1[0]);
        for (i=1; i<ndims; i++){
            sprintf(dimname, "dim1_%d", i);
            ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids1[i]);
        } 

        /* define variables */
        if (k<7){
            for (i=0; i<2; i++){
                sprintf(varname, "var0_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
            }
            for (i=2; i<nvars; i++){
                sprintf(varname, "var1_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids1, &varid[i]);
            }
        } else {
            for (i=0; i<nprocs; i++){
                sprintf(varname, "var0_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
            }
        }

        status = ncmpi_enddef(ncid);
        
	if (k == 0) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_put_vara() ...", filename);
            ncmpi_begin_indep_data(ncid);
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_put_vara(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT);
	     	TEST_HANDLE_ERR(status);
      	    }
            ncmpi_end_indep_data(ncid);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
      	} 
        if (k == 1) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_put_vara_all() ...", filename);
	    for (i=0; i<nvars; i++){
       		status = ncmpi_put_vara_all(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT);
	     	TEST_HANDLE_ERR(status);
      	    }
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
	}
	
     	if (k == 2) {
	     if (mynod == 0)
                 printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_mput_vara_all() ...", filename);
      	     status = ncmpi_mput_vara_all(ncid, nvars, varid,
                                starts_list, count_list,
                               (void **)buf, bufcount_list, datatype_list);
      	     TEST_HANDLE_ERR(status);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
     	}

        if (k == 3) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_iput_vara() and ncmpi_wait() ...", filename);
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_iput_vara(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	TEST_HANDLE_ERR(status);
                ncmpi_begin_indep_data(ncid);
	        status = ncmpi_wait(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	     	TEST_HANDLE_ERR(status);
                ncmpi_end_indep_data(ncid);
      	    }
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
      	} 
        if (k == 4) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_iput_vara() and ncmpi_wait_all() ...", filename);
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_iput_vara(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	TEST_HANDLE_ERR(status);
      	    }
	    status = ncmpi_wait_all(ncid, nvars, array_of_requests, array_of_statuses);
	    TEST_HANDLE_ERR(status);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
      	} 

        if (k == 5) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_iput_vars() and ncmpi_wait() ...", filename);
	    stride[0] = 1;
	    stride[1] = 1;
	    stride[2] = 1;
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_iput_vars(ncid, varid[i],
                            starts_list[i], count_list[i], stride,
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	TEST_HANDLE_ERR(status);
                ncmpi_begin_indep_data(ncid);
	        status = ncmpi_wait(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	     	TEST_HANDLE_ERR(status);
                ncmpi_end_indep_data(ncid);
      	    }
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
      	} 
        if (k == 6) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables into %s file by using ncmpi_iput_vars() and ncmpi_wait_all() ...", filename);
	    stride[0] = 1;
	    stride[1] = 1;
	    stride[2] = 1;
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_iput_vars(ncid, varid[i],
                            starts_list[i], count_list[i], stride,
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	TEST_HANDLE_ERR(status);
      	    }
	    status = ncmpi_wait_all(ncid, nvars, array_of_requests, array_of_statuses);
	    TEST_HANDLE_ERR(status);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
        } 
        if (k == 7) {
	    if (mynod == 0)
                printf("*** Testing to write 4 non-record variables into %s file by using ncmpi_put_var() ...", filename);
	    ncmpi_begin_indep_data(ncid);
	    status = ncmpi_put_var(ncid, varid[mynod],
	                 (const void *)(buf_var), 32, MPI_INT);
	    TEST_HANDLE_ERR(status);
	    ncmpi_end_indep_data(ncid);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
	}
        if (k == 8) {
	    if (mynod == 0)
                printf("*** Testing to write 4 non-record variables into %s file by using ncmpi_iput_var() and ncmpi_wait() ...", filename);
            i = 0;
	    status = ncmpi_iput_var(ncid, varid[mynod],
	                 (const void *)(buf_var), 32, MPI_INT, &array_of_requests[i]);

	    TEST_HANDLE_ERR(status);
	    ncmpi_begin_indep_data(ncid);
	    status = ncmpi_wait(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	    TEST_HANDLE_ERR(status);
	    ncmpi_end_indep_data(ncid);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
	}
        if (k == 9) {
	    if (mynod == 0)
                printf("*** Testing to write 4 non-record variables into %s file by using ncmpi_iput_var() and ncmpi_wait_all() ...", filename);
            i = 0;
	    status = ncmpi_iput_var(ncid, varid[mynod],
	                 (const void *)(buf_var), 32, MPI_INT, &array_of_requests[i]);

	    TEST_HANDLE_ERR(status);
	    status = ncmpi_wait_all(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	    TEST_HANDLE_ERR(status);
	    if ((mynod == 0)&&(status == NC_NOERR))
                printf("\t OK\n");                                       \
	}

    ncmpi_close(ncid);
    }/*end k_loop */

/*
    int nkeys; 
    MPI_Info_get_nkeys(info, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(info, i, key);
        MPI_Info_get_valuelen(info, key, &valuelen, &flag);
        MPI_Info_get(info, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %21s, flag = %d, valuelen = %d value = %s\n",
               i,key,flag,valuelen,value);
    }

    MPI_Info_free(&info);
*/

    MPI_Info_free(&info);
    
    for (i=0; i<nvars; i++){
        free(buf[i]);
        free(starts_list[i]);
        free(count_list[i]);
    }
    free(buf);
    free(bufcount_list);
    free(datatype_list);
    free(array_of_requests);
    free(array_of_statuses);
    free(varid);
    free(starts_list);
    free(count_list);

    MPI_Finalize();
    return 0;
}
