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
    int i, j, array_of_gsizes[3],array_of_distribs[3];
    int order, nprocs, len, **buf, bufcount, mynod;
    int array_of_dargs[3], array_of_psizes[3];
    int status;
    MPI_Offset sizes[3], array_of_starts[3];
    double write_time, *new_write_tim, write_bw;
    MPI_Offset file_size;
    double start_time, open_time, def_time, run_time;
    double *new_open_tim, *new_def_tim, *new_run_tim;
    char *pathname, filename[50];
    char dimname[20], varname[20];
    int ncid, dimids[3], rank_dim[3], *varid;
    MPI_Info info;
    MPI_Offset **starts_list, **count_list;
    int *bufcount_list;
    int ndims = 3;
    int nvars = 10;
    int k, k_loop;
    MPI_Datatype *datatype_list;
    int length;
    int mvar_flag = 0;
    NCMPI_Request *array_of_requests;
    int unlimit_flag;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynod);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    length =atoi(argv[1]);
    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    nvars = atoi(argv[2]);
    k_loop = atoi(argv[3]);
    mvar_flag = atoi(argv[4]);
    unlimit_flag = atoi(argv[5]);


/* process 0 takes the file name as a command-line argument and 
   broadcasts it to other processes */
    if (!mynod) {
	i = 1;
	while ((i < argc) && strcmp("-fname", *argv)) {
	    i++;
	    argv++;
	}
	if (i >= argc) {
	    fprintf(stderr, "\n*#  Usage: coll_perf -fname pathname\n\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	argv++;
	len = strlen(*argv);
	pathname = (char *) malloc(len+1);
        if (pathname == NULL){
		printf("pathname malloc error\n");
                return 0;		
	}	
           
	strcpy(pathname, *argv);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(pathname, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else {
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pathname = (char *) malloc(len+1);
        if (pathname == NULL){
		printf("pathname malloc error\n");
                return 0;		
	}	
	MPI_Bcast(pathname, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }


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
    bufcount_list = (int *)malloc(nvars*sizeof(int));
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
    
    array_of_requests = (NCMPI_Request *)malloc(nvars*sizeof(struct NCMPI_Req));
    
    new_open_tim = (double *)malloc(k_loop*sizeof(double));
    if (new_open_tim == NULL){
		printf("new_open_tim malloc error\n");
                return 0;		
    }	
    new_def_tim = (double *)malloc(k_loop*sizeof(double));
    if (new_def_tim == NULL){
		printf("new_def_tim malloc error\n");
                return 0;		
    }	
    new_write_tim = (double *)malloc(k_loop*sizeof(double));
    if (new_write_tim == NULL){
		printf("new_write_tim malloc error\n");
                return 0;		
    }	
    new_run_tim = (double *)malloc(k_loop*sizeof(double));
    if (new_run_tim == NULL){
		printf("new_run_tim malloc error\n");
                return 0;		
    }	
 
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


    /* mput */
    for (i=0; i<nvars; i++) {
        for (j=0; j<ndims; j++) {
           starts_list[i][j] = array_of_starts[j];
           count_list[i][j]  = length;
        }
        bufcount_list[i] = bufcount;
        datatype_list[i] = MPI_INT;
    }

    srand(mynod+1);
    for(i=0; i<nvars;i++){
        buf[i] = (int *) malloc(bufcount * sizeof(int));
        if (buf[i] == NULL){
		printf("buf[i]malloc error\n");
                return 0;		
       }	
	for (j=0; j<bufcount; j++)
		buf[i][j]=rand();
//		buf[i][j]=mynod+1;
//		buf[i][j]= mynod + 1 + 32768*i;
    }

    MPI_Info_create(&info);
    MPI_Info_set(info, "group_cyclic_fd", "enable");
//    MPI_Info_set(info, "cb_buffer_size", "1024");
    MPI_Info_set(info, "cb_buffer_size", "16777216");
/*    MPI_Info_set(info, "romio_no_indep_rw", "true");*/
    MPI_Info_set(info, "romio_cb_write", "true");

    for (k=0; k<k_loop; k++){
      sprintf(filename, "%s.%d.%d.%d.%d.nc", pathname, length, nvars, mvar_flag, k);
      MPI_Barrier(MPI_COMM_WORLD);
      start_time = MPI_Wtime();
      status = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_OFFSET,
                        info, &ncid);
      TEST_HANDLE_ERR(status);
      /* define dimensions */
      if (unlimit_flag == 1) {
        sprintf(dimname, "dim_%d", 0);
        ncmpi_def_dim(ncid, dimname, NC_UNLIMITED, &dimids[0]);
        for (i=1; i<ndims; i++){
         sprintf(dimname, "dim_%d", i);
         ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids[i]);
       } 
      } else {
        for (i=0; i<ndims; i++){
         sprintf(dimname, "dim_%d", i);
         ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids[i]);
       } 
      }

      /* define variables */
      for (i=0; i<nvars; i++){
        sprintf(varname, "var_%d", i);
        ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids, &varid[i]);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      open_time = MPI_Wtime()-start_time;
    
      status = ncmpi_enddef(ncid);
      MPI_Barrier(MPI_COMM_WORLD);
      def_time = MPI_Wtime()-start_time-open_time;

/* to eliminate paging effects, do the operations once but don't time
   them */
        if (mvar_flag == 0) {
	      for (i=0; i<nvars; i++){
       		  status = ncmpi_put_vara_all(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT);
	     	  TEST_HANDLE_ERR(status);
      	      }
      	} 
        if (mvar_flag == 1) {
//      	     status = ncmpi_put_mvara_all_record(ncid, nvars, varid,
      	     status = ncmpi_put_mvara_all(ncid, nvars, varid,
                                starts_list, count_list,
                               (const void **)buf, bufcount_list, datatype_list);
      	     TEST_HANDLE_ERR(status);
     	}
        if (mvar_flag == 2) {
	      for (i=0; i<nvars; i++){
       		  status = ncmpi_iput_vara_all(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	  TEST_HANDLE_ERR(status);
	          ncmpi_wait(array_of_requests[i]);
      	      }
      	} 
        if (mvar_flag == 3) {
	      for (i=0; i<nvars; i++){
       		  status = ncmpi_iput_vara_all(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	  TEST_HANDLE_ERR(status);
      	      }
	      ncmpi_waitall(nvars, array_of_requests);
      	} 
 
 
      MPI_Barrier(MPI_COMM_WORLD);
      write_time = MPI_Wtime() - start_time - open_time - def_time;
//      printf("mynod:%d, write_time:%.6f\n", mynod, write_time);
/*    ncmpi_get_file_info(ncid, &info); */
    
      ncmpi_close(ncid);
      run_time = MPI_Wtime() - start_time;
      MPI_Allreduce(&open_time, &new_open_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&def_time, &new_def_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&write_time, &new_write_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&run_time, &new_run_tim[k], 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
    }/*end k_loop */

    if (mynod == 0) {
      int top = 0;
      for (k=0; k<(k_loop-1); k++){
        if (new_run_tim[top]>new_run_tim[k+1]) top = k+1;
      }

      file_size = 1024+nvars*array_of_gsizes[0]*array_of_gsizes[1]*array_of_gsizes[2]*sizeof(int);
      write_bw = file_size/new_run_tim[top]/1024.0/1024.0;
      fprintf(stderr, "mvar nvars:%d, Global array size %d x %d x %d integers, local array size: %d x %d x%d, filesize:%d\n", nvars, array_of_gsizes[0], array_of_gsizes[1], array_of_gsizes[2],sizes[0], sizes[1], sizes[2], file_size);
      fprintf(stderr, "%dx%dx%d, %d: nvars:%d, loop:%d, k:%d, open_t = %f, def_t =%f, write_t = %f sec,run_t = %f sec, filesize: %ld, bandwidth = %f Mbytes/sec\n", sizes[0], sizes[1], sizes[2],mvar_flag, nvars, k_loop, top, new_open_tim[top], new_def_tim[top], new_write_tim[top], new_run_tim[top],file_size, write_bw); 
    }
    
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

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Info_free(&info);
    
    free(new_open_tim);
    free(new_def_tim);
    free(new_write_tim);
    free(new_run_tim);

    for (i=0; i<nvars; i++){
    free(buf[i]);
    free(starts_list[i]);
    free(count_list[i]);
    }
    free(buf);
    free(bufcount_list);
    free(varid);
    free(starts_list);
    free(count_list);
    free(pathname);

    MPI_Finalize();
    return 0;
}
