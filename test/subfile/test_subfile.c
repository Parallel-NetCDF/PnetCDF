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

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    extern char *optarg;
    extern int optind;
    char dimname[20], varname[20];
    int i, j, err, nerrs=0, verbose=0, nprocs, rank, ncid, *varid=NULL;
    int ndims=3, ngatts, unlimdimid, dimids0[3], rank_dim[3];
    int num_files, **buf, array_of_psizes[3], array_of_gsizes[3];
    MPI_Offset bufcount, array_of_starts[3], **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    MPI_Datatype *datatype_list;
    MPI_Info info_used=MPI_INFO_NULL;
    double stim, write_tim, new_write_tim, write_bw, read_bw;
    double open_tim, new_open_tim, read_tim, new_read_tim;
    double close_tim, new_close_tim;

    int num_sf = 2;
    int par_dim_id = 0; /* default is 0 */
    int do_read = 0;
    int nvars = 1;
    int length = 8;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    buf = (int **)malloc(sizeof(int*) * nvars);
    if (buf == NULL){
        printf("buf malloc error\n");
        nerrs++; goto fn_exit;
    }
    bufcount_list = (MPI_Offset *)malloc(sizeof(MPI_Offset)*nvars);
    if (bufcount_list == NULL){
        printf("bufcount_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    starts_list = (MPI_Offset **)malloc(sizeof(MPI_Offset *)*nvars);
    if (starts_list== NULL){
        printf("starts_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    count_list = (MPI_Offset **)malloc(sizeof(MPI_Offset *)*nvars);
    if (count_list == NULL){
        printf("count_list malloc error\n");
        nerrs++; goto fn_exit;
    }
    datatype_list = (MPI_Datatype*)malloc(sizeof(MPI_Datatype)*nvars);
    if (datatype_list == NULL){
        printf("count_list malloc error\n");
        nerrs++; goto fn_exit;
    }

    for (i=0; i<nvars; i++) {
        starts_list[i] = (MPI_Offset *)malloc(sizeof(MPI_Offset)*ndims);
        if (starts_list[i] == NULL){
            printf("starts_list[%d] malloc error\n", i);
            nerrs++; goto fn_exit;
        }
        count_list[i] = (MPI_Offset *)malloc(sizeof(MPI_Offset)*ndims);
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
        buf[i] = (int *) malloc(sizeof(int) * bufcount);
        if (buf[i] == NULL){
            printf("buf[i]malloc error\n");
            nerrs++; goto fn_exit;
        }

        for (j=0; j<bufcount; j++)
            buf[i][j]=rank+1;
    }

    /* set all non-record variable to be subfiled */
    char tmp[10];
    sprintf(tmp, "%d", num_sf);
    MPI_Info_set(info, "nc_num_subfiles", tmp);
    MPI_Info_set(info, "pnetcdf_subfiling", "enable");

    if (do_read == 1) goto read;

    stim = MPI_Wtime();

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
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
    varid = (int *)malloc(sizeof(int)*nvars);
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
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);
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

    int nfiles, ncids[10];

    /* NULL argument test */
    err = ncmpi_inq_files_opened(NULL, NULL);
    EXP_ERR(NC_EINVAL)

    /* check if there are files still left opened */
    err = ncmpi_inq_files_opened(&nfiles, ncids);
    CHECK_ERR
    if (nfiles > 0) printf("nfiles %d still opened\n",nfiles);

    /* open the subfiles to validate the file format */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank > 0) goto fn_exit;

    num_files = (nprocs < num_sf) ? 1 : num_sf;
    for (i=0; i<num_files; i++) {
        char filename[512];
        sprintf(filename, "%s.subfile_%d.nc", out_path, i);
        err = ncmpi_open(MPI_COMM_SELF, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
        CHECK_ERR
        err = ncmpi_close(ncid);
        CHECK_ERR
        unlink(filename);
    }

fn_exit:
    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "subfiling", opt, test_io);

    MPI_Finalize();

    return err;
}
