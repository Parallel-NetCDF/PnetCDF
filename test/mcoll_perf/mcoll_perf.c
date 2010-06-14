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

#define HANDLE_ERROR {                                \
    if (status != NC_NOERR)                           \
        printf("Error at line %d (%s)\n", __LINE__,   \
               ncmpi_strerror(status));               \
}

#define HANDLE_DIFF(str) {                                       \
    int doStop, isDiff = (str[0] == '\0') ? 0 : 1;               \
    MPI_Allreduce(&isDiff, &doStop, 1, MPI_INT, MPI_MAX, comm);  \
    if (doStop) {                                                \
        printf("P%d: diff at %s",rank,str);                      \
        MPI_Finalize();                                          \
        return 1;                                                \
    }                                                            \
}

#define CHECK_GLOBAL_ATT_DIFF(type, func, nctype) {                          \
    int   pos, len = attlen1 * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, NC_GLOBAL, name1, b1);                              \
    HANDLE_ERROR                                                             \
    status = func(ncid2, NC_GLOBAL, name2, b2);                              \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"attribute[%d] %s: %s buf1 != buf2 at position %d\n",    \
                i,name1,#nctype,pos);                                        \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_ATT_DIFF(type, func, nctype) {                             \
    int   pos, len = attlen1 * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, i, name1, b1);                                      \
    HANDLE_ERROR                                                             \
    status = func(ncid2, i, name2, b2);                                      \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"variable[%d] %s: attribute[%d] %s: %s buf1 != buf2 at position %d\n", \
                i,name,j,name1,#nctype,pos);                                 \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}


#define CHECK_VAR_DIFF(type, func, nctype) {                                 \
    int   pos, len = varsize * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, i, start, shape, b1);                               \
    HANDLE_ERROR                                                             \
    status = func(ncid2, i, start, shape, b2);                               \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"variable[%d] %s: %s buf1 != buf2 at position %d\n",     \
                i,name,#nctype,pos);                                        \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}


int ncmpi_diff(char *filename1, char *filename2) {
    int i, j, status, isRecvar, rank, nprocs;
    int ncid1, ndims1, nvars1, natts1, unlimdimid1, dimids1[NC_MAX_DIMS];
    int ncid2, ndims2, nvars2, natts2, unlimdimid2, dimids2[NC_MAX_DIMS];
    char str[512], name1[NC_MAX_NAME], name2[NC_MAX_NAME], name[NC_MAX_NAME];
    MPI_Offset shape[NC_MAX_VAR_DIMS], varsize, start[NC_MAX_VAR_DIMS];
    MPI_Offset attlen1, dimlen1, attlen2, dimlen2;
    ncmpi_type type1, type2;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    str[0] = '\0';

    status = ncmpi_open(comm, filename1, NC_NOWRITE, MPI_INFO_NULL, &ncid1);
    HANDLE_ERROR
    status = ncmpi_open(comm, filename2, NC_NOWRITE, MPI_INFO_NULL, &ncid2);
    HANDLE_ERROR

    /**
     * Inquire the dataset definitions of input dataset AND
     * Add dataset definitions for output dataset.
     */
    status = ncmpi_inq(ncid1, &ndims1, &nvars1, &natts1, &unlimdimid1);
    HANDLE_ERROR
    status = ncmpi_inq(ncid2, &ndims2, &nvars2, &natts2, &unlimdimid2);
    HANDLE_ERROR
    if (ndims1 != ndims2)
        sprintf(str,"ndims1(%d) != ndims2(%d)\n",ndims1, ndims2);
    HANDLE_DIFF(str)
    if (nvars1 != nvars2)
        sprintf(str,"nvars1(%d) != nvars2(%d)\n",nvars1, nvars2);
    HANDLE_DIFF(str)
    if (natts1 != natts2)
        sprintf(str,"natts1(%d) != natts2(%d)\n",natts1, natts2);
    HANDLE_DIFF(str)

    /* Inquire global attributes, assume CHAR attributes. */
    for (i=0; i<natts1; i++) {
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name1);
        HANDLE_ERROR
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0)
            sprintf(str,"attribute[%d] name1(%s) != name2(%s)\n",i,name1,name2);
        HANDLE_DIFF(str)

        status = ncmpi_inq_att(ncid1, NC_GLOBAL, name1, &type1, &attlen1);
        HANDLE_ERROR
        status = ncmpi_inq_att(ncid2, NC_GLOBAL, name2, &type2, &attlen2);
        HANDLE_ERROR
        if (type1 != type2)
            sprintf(str,"attribute[%d] %s: type1(%d) != type2(%d)\n",i,name1,type1,type2);
        HANDLE_DIFF(str)
        if (attlen1 != attlen2)
            sprintf(str,"attribute[%d] %s: attlen1(%lld) != attlen2(%lld)\n",i,name1,(long long int) attlen1,(long long int)attlen2);
        HANDLE_DIFF(str)
        switch (type1) {
            case NC_CHAR:   CHECK_GLOBAL_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
            case NC_SHORT:  CHECK_GLOBAL_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
            case NC_INT:    CHECK_GLOBAL_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
            case NC_FLOAT:  CHECK_GLOBAL_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
            case NC_DOUBLE: CHECK_GLOBAL_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
            default: ; /* TODO: handle unexpected types */
        }
    }

    /* Inquire dimension */
    for (i=0; i<ndims1; i++) {
        status = ncmpi_inq_dim(ncid1, i, name1, &dimlen1);
        HANDLE_ERROR
        status = ncmpi_inq_dim(ncid2, i, name2, &dimlen2);
        HANDLE_ERROR
        if (dimlen1 != dimlen2)
            sprintf(str,"dimension[%d] %s: dimlen1(%lld) != dimlen2(%lld)\n",i,name1,(long long int)dimlen1,(long long int)dimlen2);
        HANDLE_DIFF(str)
    }

    /* Inquire variables */
    for (i=0; i<nvars1; i++) {
        status = ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
        HANDLE_ERROR
        status = ncmpi_inq_var(ncid2, i, name2, &type2, &ndims2, dimids2, &natts2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0)
            sprintf(str,"variable[%d]: name1(%s) != name2(%s)\n",i,name1,name2);
        HANDLE_DIFF(str)
        if (type1 != type2)
            sprintf(str,"variable[%d] %s: type1(%d) != type2(%d)\n",i,name1,type1,type2);
        HANDLE_DIFF(str)
        if (ndims1 != ndims2)
            sprintf(str,"variable[%d] %s: ndims1(%d) != ndims2(%d)\n",i,name1,ndims1,ndims2);
        HANDLE_DIFF(str)
        for (j=0; j<ndims1; j++) {
            if (dimids1[j] != dimids2[j])
                sprintf(str,"variable[%d] %s: dimids1[%d]=%d != dimids2[%d]=%d\n",i,name1,j,dimids1[j],j,dimids2[j]);
            HANDLE_DIFF(str)
        }
        if (natts1 != natts2)
            sprintf(str,"variable[%d] %s: natts1(%d) != natts2(%d)\n",i,name1,natts1,natts2);
        HANDLE_DIFF(str)

        strcpy(name,name1);

        /* var attributes, assume CHAR attributes */
        for (j=0; j<natts1; j++) {
            status = ncmpi_inq_attname(ncid1, i, j, name1);
            HANDLE_ERROR
            status = ncmpi_inq_attname(ncid2, i, j, name2);
            HANDLE_ERROR
            if (strcmp(name1, name2) != 0)
                sprintf(str,"variable[%d] %s: attr name[%d] (%s) != (%s)\n",i,name,j,name1,name2);
            HANDLE_DIFF(str)

            status = ncmpi_inq_att(ncid1, i, name1, &type1, &attlen1);
            HANDLE_ERROR
            status = ncmpi_inq_att(ncid2, i, name2, &type2, &attlen2);
            HANDLE_ERROR
            if (type1 != type2)
                sprintf(str,"variable[%d] %s: attr type[%d] (%d) != (%d)\n",i,name,j,type1,type2);
            HANDLE_DIFF(str)
            if (attlen1 != attlen2)
                sprintf(str,"variable[%d] %s: attr attlen[%d] (%lld) != (%lld)\n",i,name,j,(long long int)attlen1,(long long int)attlen2);
            HANDLE_DIFF(str)

            switch (type1) {
                case NC_CHAR:   CHECK_VAR_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
                case NC_SHORT:  CHECK_VAR_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
                case NC_INT:    CHECK_VAR_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
                case NC_FLOAT:  CHECK_VAR_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
                case NC_DOUBLE: CHECK_VAR_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
                default: ; /* TODO: handle unexpected types */
            }
        }
    }

    /**
     * Read data of variables from input dataset 
     * (ONLY DEAL WITH: NC_INT, NC_FLOAT, NC_DOUBLE for now)
     * Write the data out to the corresponding variables in the output dataset
     *
     *  Data Partition (Assume 4 processors):
     *   square: 2-D, (Block, *), 25*100 from 100*100
     *   cube:   3-D, (Block, *, *), 25*100*100 from 100*100*100
     *   xytime: 3-D, (Block, *, *), 25*100*100 from 100*100*100
     *   time:   1-D, Block-wise, 25 from 100
     *
     *  Data Mode API: collective
     */

    for (i=0; i<NC_MAX_VAR_DIMS; i++)
        start[i] = 0;

    for (i=0; i<nvars1; i++) {
        isRecvar = 0;
        varsize = 1;
        ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
        strcpy(name,name1);
        for (j=0; j<ndims1; j++) {
            status = ncmpi_inq_dim(ncid1, dimids1[j], name2, shape + j);
            HANDLE_ERROR
            /* name2 will be discarded */
            if (j == 0) {
                shape[j] /= nprocs;
                start[j] = shape[j] * rank;
            }
            varsize *= shape[j];
            if (dimids1[j] == unlimdimid1)
                isRecvar = 1;
        }
        switch (type1) {
            case NC_CHAR:   CHECK_VAR_DIFF(char,   ncmpi_get_vara_text_all,   NC_CHAR)
            case NC_SHORT:  CHECK_VAR_DIFF(short,  ncmpi_get_vara_short_all,  NC_SHORT)
            case NC_INT:    CHECK_VAR_DIFF(int,    ncmpi_get_vara_int_all,    NC_INT)
            case NC_FLOAT:  CHECK_VAR_DIFF(float,  ncmpi_get_vara_float_all,  NC_FLOAT)
            case NC_DOUBLE: CHECK_VAR_DIFF(double, ncmpi_get_vara_double_all, NC_DOUBLE)
            default: ; /* TODO: handle unexpected types */
        }
    }

    status = ncmpi_close(ncid1);
    HANDLE_ERROR
    status = ncmpi_close(ncid2);
    HANDLE_ERROR

    return NC_NOERR;
}


    
int main(int argc, char **argv)
{
    int i, j, array_of_gsizes[3],array_of_distribs[3];
    int order, nprocs, len, **buf, mynod;
    MPI_Offset bufcount;
    int array_of_dargs[3], array_of_psizes[3];
    int status;
    MPI_Offset sizes[3], array_of_starts[3], stride[3];
    char basename[50], filename[100];
    char filename1[100], filename2[100], filename3[100];
    char dimname[20], varname[20];
    int ncid, dimids0[3], dimids1[3], rank_dim[3], *varid;
    MPI_Info info;
    MPI_Offset **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    int ndims = 3;
    int nvars = 10;
    int k;
    MPI_Datatype *datatype_list;
    int length;
    int *array_of_requests;
    int *array_of_statuses;
    int *buf_var;
//    int buf_var[32] ={1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2,
//                      3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4};

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynod);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    length = 2;
    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    nvars = 4;
//    strcpy(basename,  "/pvfs2/kgao/test_non_blocking/test");
    strcpy(basename,  "./test");

/* process 0 takes the file name as a command-line argument and
   broadcasts it to other processes */
   if (argc > 1){
     if (!mynod) {
        i = 1;
        while ((i < argc) && strcmp("-fname", *argv)) {
            i++;
            argv++;
        }
        if (i >= argc) {
            fprintf(stderr, "\n*#  Usage: mcoll_perf -fname pathname\n\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        argv++;
        len = strlen(*argv);
//        pathname = (char *) xmalloc(len+1);

        strcpy(basename, *argv);
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(basename, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      }
      else {
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        pathname = (char *) xmalloc(len+1);
        MPI_Bcast(basename, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      }
    }
    order = MPI_ORDER_C;
    
    buf = (int **)malloc(nvars*sizeof(int*));
    if (buf == NULL){
	printf("buf malloc error\n");
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
    buf_var = (int *) malloc(bufcount*nprocs*sizeof(int));
    for (i=0; i<bufcount*nprocs; i++)
        buf_var[i] = mynod + 1;  

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
	if (k==0)
	    strcpy(filename1, filename);
        else if (k==7)
	    strcpy(filename2, filename);
	else 
	    strcpy(filename3, filename);

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
            varid = (int *)malloc(nvars*sizeof(int));
            if (varid == NULL){
	        printf("varid malloc error\n");
                return 0;		
            }	
            for (i=0; i<2; i++){
                sprintf(varname, "var0_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
            }
            for (i=2; i<nvars; i++){
                sprintf(varname, "var1_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids1, &varid[i]);
            }
        } else {
            varid = (int *)malloc(nprocs*sizeof(int));
            if (varid == NULL){
	        printf("varid malloc error\n");
                return 0;		
            }	
            for (i=0; i<nprocs; i++){
                sprintf(varname, "var0_%d", i);
                ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
            }
        }

        status = ncmpi_enddef(ncid);
        
        if (k == 0) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_put_vara_all() ...");
	    for (i=0; i<nvars; i++){
       		status = ncmpi_put_vara_all(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT);
	     	TEST_HANDLE_ERR(status);
      	    }
	}
	if (k == 1) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_put_vara() ...");
            ncmpi_begin_indep_data(ncid);
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_put_vara(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT);
	     	TEST_HANDLE_ERR(status);
      	    }
            ncmpi_end_indep_data(ncid);
      	} 
     	if (k == 2) {
	     if (mynod == 0)
                 printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_mput_vara_all() ...");
      	     status = ncmpi_mput_vara_all(ncid, nvars, varid,
                                starts_list, count_list,
                               (void **)buf, bufcount_list, datatype_list);
      	     TEST_HANDLE_ERR(status);
     	}

        if (k == 3) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_iput_vara() and ncmpi_wait() ...");
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
      	} 
        if (k == 4) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_iput_vara() and ncmpi_wait_all() ...");
	    for (i=0; i<nvars; i++){
       	        status = ncmpi_iput_vara(ncid, varid[i],
                            starts_list[i], count_list[i],
                            (const void *)&(buf[i][0]), bufcount_list[i], MPI_INT, &array_of_requests[i]);
	     	TEST_HANDLE_ERR(status);
      	    }
	    status = ncmpi_wait_all(ncid, nvars, array_of_requests, array_of_statuses);
	    TEST_HANDLE_ERR(status);
      	} 

        if (k == 5) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_iput_vars() and ncmpi_wait() ...");
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
      	} 
        if (k == 6) {
	    if (mynod == 0)
                printf("*** Testing to write 2 non-record variables and 2 record variables by using ncmpi_iput_vars() and ncmpi_wait_all() ...");
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
        } 
        if (k == 7) {
	    if (mynod == 0)
                printf("*** Testing to write %d non-record variable(s) by using ncmpi_put_var() ...", nprocs);
	    ncmpi_begin_indep_data(ncid);
	    status = ncmpi_put_var(ncid, varid[mynod],
	                 (const void *)(buf_var), bufcount*nprocs, MPI_INT);
	    TEST_HANDLE_ERR(status);
	    ncmpi_end_indep_data(ncid);
	}
        if (k == 8) {
	    if (mynod == 0)
                printf("*** Testing to write %d non-record variable(s) by using ncmpi_iput_var() and ncmpi_wait() ...", nprocs);
            i = 0;
	    status = ncmpi_iput_var(ncid, varid[mynod],
	                 (const void *)(buf_var), bufcount*nprocs, MPI_INT, &array_of_requests[i]);

	    TEST_HANDLE_ERR(status);
	    ncmpi_begin_indep_data(ncid);
	    status = ncmpi_wait(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	    TEST_HANDLE_ERR(status);
	    ncmpi_end_indep_data(ncid);
	}
        if (k == 9) {
	    if (mynod == 0)
                printf("*** Testing to write %d non-record variable(s) by using ncmpi_iput_var() and ncmpi_wait_all() ...", nprocs);
            i = 0;
	    status = ncmpi_iput_var(ncid, varid[mynod],
	                 (const void *)(buf_var), bufcount*nprocs, MPI_INT, &array_of_requests[i]);

	    TEST_HANDLE_ERR(status);
	    status = ncmpi_wait_all(ncid, 1, &array_of_requests[i], &array_of_statuses[i]);
	    TEST_HANDLE_ERR(status);
	}
        ncmpi_close(ncid);
        if (status == NC_NOERR){
            if ((k>0)&&(k<7)){
	        status = ncmpi_diff(filename1, filename3);
	        if ((mynod == 0)&&(status == NC_NOERR))
                    printf("\t OK\n");                                       
             } else if (k>7){
	         status = ncmpi_diff(filename2, filename3);
	         if ((mynod == 0)&&(status == NC_NOERR))
                     printf("\t OK\n");                                       
	     } else {
	         if (mynod == 0)
                 printf("\t OK\n");                                       
             }
        }
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

    MPI_Info_free(&info);
    
    for (i=0; i<nvars; i++){
        free(buf[i]);
        free(starts_list[i]);
        free(count_list[i]);
    }
    free(buf);
    free(buf_var);
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
