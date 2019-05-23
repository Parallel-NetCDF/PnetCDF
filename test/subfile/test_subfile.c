/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

/* The file name is taken as a command-line argument. */

/* Measures the I/O bandwidth for writing/reading a 3D
   block-distributed array to a file corresponding to the global array
   in row-major (C) order.
   Note that the file access pattern is noncontiguous.

   Array size 128^3. For other array sizes, change array_of_gsizes below.
*/

int main(int argc, char **argv)
{
    int opt, verbose=0;
    extern char *optarg;
    extern int optind;
    int i, j, array_of_gsizes[3];
    int nprocs, len, **buf, rank;
    MPI_Offset bufcount;
    int array_of_psizes[3];
    int err;
    MPI_Offset array_of_starts[3];
    char *fbasename=NULL;
    char dimname[20], varname[20];
    int ncid, dimids0[3], rank_dim[3], *varid=NULL;
    MPI_Info info=MPI_INFO_NULL, info_used=MPI_INFO_NULL;
    MPI_Offset **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    int ndims=3, nvars=1, ngatts, unlimdimid;
    MPI_Datatype *datatype_list;
    int length = 8;
    double stim, write_tim, new_write_tim, write_bw;
    double read_tim, new_read_tim, read_bw;
    double open_tim, new_open_tim;
    double close_tim, new_close_tim;
    int num_sf = 2;
    int par_dim_id = 0; /* default is 0 */
    int do_read = 0;
    int nerrs=0;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* process 0 takes the file name as a command-line argument and
       broadcasts it to other processes */
    if (rank == 0) {
	while ((opt = getopt(argc, argv, "f:s:p:n:l:r")) != EOF) {
	    switch (opt) {
	    case 'f': fbasename = optarg;
		break;
	    case 's': num_sf = (int)strtol(optarg,NULL,10);
		break;
	    case 'r': do_read = 1;
		break;
            case 'p': par_dim_id = (int)strtol(optarg,NULL,10);
                break;
            case 'n': nvars = (int)strtol(optarg,NULL,10);
                break;
            case 'l': length = (int)strtol(optarg,NULL,10);
                break;
	    default:
		break;
	    }
	}
	if (fbasename == NULL) {
	    fprintf(stderr, "\n*#  Usage: test_subfile -f pathname -s num_sf -p par_dim_id \n\n");
	    nerrs++;
	}
    }
    MPI_Bcast(&nerrs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (nerrs > 0) {
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
	len = (fbasename == NULL) ?  0 : strlen(fbasename);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	fbasename = (char *) malloc(len+1);
    }
    MPI_Bcast(fbasename, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_sf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&par_dim_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nvars, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&do_read, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for subfiling", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    buf = (int **)malloc(nvars*sizeof(int*));
    if (buf == NULL){
        printf("buf malloc error\n");
        nerrs++; goto fn_exit;
    }
    bufcount_list = (MPI_Offset *)malloc(nvars*sizeof(MPI_Offset));
    if (bufcount_list == NULL){
        printf("bufcount_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    starts_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (starts_list== NULL){
        printf("starts_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    count_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (count_list == NULL){
        printf("count_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    datatype_list = (MPI_Datatype*)malloc(nvars*sizeof(MPI_Datatype));
    if (datatype_list == NULL){
        printf("count_list malloc error\n");
        nerrs++; goto fn_exit;
    }

    for (i=0; i<nvars; i++) {
	starts_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	if (starts_list[i] == NULL){
	    printf("starts_list[%d] malloc error\n", i);
	    nerrs++; goto fn_exit;
	}
	count_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	if (count_list[i] == NULL){
	    printf("count_list[%d] malloc error\n", i);
	    nerrs++; goto fn_exit;
	}
    }

    bufcount = 1;
    for (i=0; i<ndims; i++) {
        array_of_psizes[i] = 0;
        bufcount *= length;
    }
    MPI_Dims_create(nprocs, ndims, array_of_psizes);
    if (verbose && rank == 0)
        for(i=0; i<ndims; i++)
	    printf("array_of_psizes[%d]=%d\n", i, array_of_psizes[i]);

    /* subarray in each process is len x len x len */
    for (i=0; i<ndims; i++)
        array_of_gsizes[i] = length * array_of_psizes[i];

    /* mynd's process rank in each dimension (in MPI_ORDER_C) */
    rank_dim[2] =  rank %  array_of_psizes[2];
    rank_dim[1] = (rank /  array_of_psizes[2]) % array_of_psizes[1];
    rank_dim[0] =  rank / (array_of_psizes[2]  * array_of_psizes[1]);

    /* starting coordinates of the subarray in each dimension */
    for (i=0; i<ndims; i++)
        array_of_starts[i] = (MPI_Offset)length * rank_dim[i];

    for (i=0; i<nvars; i++) {
	for (j=0; j<ndims; j++) {
	    starts_list[i][j] = array_of_starts[j];
	    count_list[i][j]  = length;
	}
        bufcount_list[i] = bufcount;
        datatype_list[i] = MPI_INT;
    }

    for (i=0; i<nvars; i++) {
	buf[i] = (int *) malloc(bufcount * sizeof(int));
	if (buf[i] == NULL){
	    printf("buf[i]malloc error\n");
	    nerrs++; goto fn_exit;
	}

	for (j=0; j<bufcount; j++)
	    buf[i][j]=rank+1;
    }

    MPI_Info_create(&info);
    /* set all non-record variable to be subfiled */
    char tmp[10];
    sprintf(tmp, "%d", num_sf);
    MPI_Info_set(info, "nc_num_subfiles", tmp);
    MPI_Info_set(info, "pnetcdf_subfiling", "enable");

    if (do_read == 1) goto read;

    stim = MPI_Wtime();
    err = ncmpi_create(MPI_COMM_WORLD, fbasename, NC_CLOBBER|NC_64BIT_DATA,
                       info, &ncid);
    CHECK_ERR

    open_tim = MPI_Wtime() - stim;

    MPI_Allreduce(&open_tim, &new_open_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    if (verbose && rank == 0)
        printf("create time = %f sec\n", new_open_tim);

    /* define dimensions */
    for (i=0; i<ndims; i++){
        sprintf(dimname, "dim0_%d", i);
        err = ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids0[i]);
        CHECK_ERR
    }

    /* define variables */
    varid = (int *)malloc(nvars*sizeof(int));
    for (i=0; i<nvars; i++) {
	sprintf(varname, "var0_%d", i);
	err = ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
	CHECK_ERR
    }

    if (par_dim_id != 0) {
        for (i=0; i<nvars; i++) {
            err = ncmpi_put_att_int(ncid, varid[i], "par_dim_id",
                                    NC_INT, 1, &dimids0[par_dim_id]);
	    CHECK_ERR
        }
    }

    /* set all non-record variable to be subfiled */
    /*
    MPI_Info_set(info, "nc_num_subfiles", "2");
    err = ncmpi_set_var_info(ncid, varid, info);
    CHECK_ERR
    */

    err = ncmpi_enddef(ncid);
    CHECK_ERR

    /* test ncmpi_inq_var() */
    for (i=0; i<nvars; i++) {
        char name[128];
        nc_type typep;
        int ndimsp, dimids[3], nattsp;

        err = ncmpi_inq_var(ncid, varid[i], name, &typep, &ndimsp, dimids,
                            &nattsp);
        CHECK_ERR

	sprintf(varname, "var0_%d", i);
        if (strcmp(name, varname)) {
            printf("Error at line %d in %s: unexpected var[%d] name %s, should be %s\n",
            __LINE__,__FILE__,i,name,varname);
            nerrs++;
            continue;
        }
        if (typep != NC_INT) {
            printf("Error at line %d in %s: unexpected var[%d] type %d, should be %d\n",
            __LINE__,__FILE__,i,typep,NC_INT);
            nerrs++;
            continue;
        }
        if (ndimsp != ndims) {
            printf("Error at line %d in %s: unexpected var[%d] ndims %d, should be %d\n",
            __LINE__,__FILE__,i,ndimsp,ndims);
            nerrs++;
            continue;
        }
        for (j=0; j<ndims; j++) {
            if (dimids[j] != dimids0[j]) {
                printf("Error at line %d in %s: unexpected var[%d] dimids[%d] %d, should be %d\n",
            __LINE__,__FILE__,i,j,dimids0[j],dimids[j]);
                nerrs++;
                continue;
            }
        }
        /*
        printf("var[%d] %s has %d attributes\n",i,name,nattsp);
        */
    }

#if 0
    if (rank == 0)
        printf("*** Testing to write 1 non-record variable by using ncmpi_put_vara_all() ...");
#endif
    stim = MPI_Wtime();
    for (i=0; i<nvars; i++) {
        err = ncmpi_put_vara_all(ncid, varid[i],
                                 starts_list[i], count_list[i],
                                 buf[i],
                                 bufcount_list[i], MPI_INT);
        CHECK_ERR
    }
    write_tim = MPI_Wtime() - stim;

    MPI_Allreduce(&write_tim, &new_write_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (verbose && rank == 0) {
        write_bw = ((double)array_of_gsizes[0]*(double)array_of_gsizes[1]*(double)array_of_gsizes[2]*(double)sizeof(int)*(double)nvars)/(new_write_tim*1024.0*1024.0);
        printf("Global array size %d x %d x %d integers\n", array_of_gsizes[0], array_of_gsizes[1], array_of_gsizes[2]);
        printf("Collective write time = %f sec, Collective write bandwidth = %f Mbytes/sec\n", new_write_tim, write_bw);
    }

    err = ncmpi_inq_file_info(ncid, &info_used);
    CHECK_ERR

    stim = MPI_Wtime();
    err = ncmpi_close(ncid);
    CHECK_ERR
    close_tim = MPI_Wtime() - stim;

    MPI_Allreduce(&close_tim, &new_close_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (verbose && rank == 0) {
        fprintf(stderr, "close time = %f sec\n", new_close_tim);
    }

    goto end;

read:
    err = ncmpi_open(MPI_COMM_WORLD, fbasename, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    stim = MPI_Wtime();

    /**
     * Inquire the dataset definitions of input dataset AND
     * Add dataset definitions for output dataset.
     */

    err = ncmpi_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
    CHECK_ERR

    for (i=0; i<nvars; i++) {
        err = ncmpi_get_vara_all(ncid, i,
                                 starts_list[i], count_list[i],
                                 buf[i], bufcount_list[i], MPI_INT);
        CHECK_ERR
    }
    read_tim = MPI_Wtime() - stim;

    MPI_Allreduce(&read_tim, &new_read_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (verbose && rank == 0) {
        read_bw = ((double)array_of_gsizes[0]*(double)array_of_gsizes[1]*(double)array_of_gsizes[2]*sizeof(int)*(double)nvars)/(new_read_tim*1024.0*1024.0);
        printf("Collective read time = %f sec, Collective read bandwidth = %f Mbytes/sec\n", new_read_tim, read_bw);
    }

    err = ncmpi_inq_file_info(ncid, &info_used);
    CHECK_ERR

    err = ncmpi_close(ncid);
    CHECK_ERR

end:
    if (info      != MPI_INFO_NULL) MPI_Info_free(&info);
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);

    for (i=0; i<nvars; i++){
        free(buf[i]);
        free(starts_list[i]);
        free(count_list[i]);
    }
    free(buf);
    free(bufcount_list);
    free(datatype_list);
    if (!do_read) free(varid);
    free(starts_list);
    free(count_list);
    if (rank > 0) free(fbasename);

    MPI_Offset malloc_size, sum_size;
    int nfiles, ncids[10];

    /* NULL argument test */
    err = ncmpi_inq_files_opened(NULL, NULL);
    EXP_ERR(NC_EINVAL)

    /* check if there are files still left opened */
    err = ncmpi_inq_files_opened(&nfiles, ncids);
    CHECK_ERR
    if (nfiles > 0) printf("nfiles %d still opened\n",nfiles);

    /* check for any PnetCDF internal malloc residues */
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}
