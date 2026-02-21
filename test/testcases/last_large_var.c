/*********************************************************************
 *
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests the special case when there is no record variable, the
 * last fixed-size variable can be larger than 2GiB if its starting file offset
 * is less than 2GiB. See NetCDF Classic Format Limitations (The NetCDF Users
 * Guide).  Quoted here:
 * "If you don't use the unlimited dimension, only one variable can exceed 2
 * GiB in size, but it can be as large as the underlying file system permits.
 * It must be the last variable in the dataset, and the offset to the beginning
 * of this variable must be less than about 2 GiB."
 * http://www.unidata.ucar.edu/software/netcdf/old_docs/docs_3_6_3/netcdf-c/nc_005fcreate.html
 *
 *    To compile:
 *        mpicc -O2 last_large_var.c -o last_large_var -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * NC file produced by this example program:
 *
 *    % mpiexec -n 4 ./last_large_var /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    netcdf testfile {
 *    dimensions:
 *    	    Y = 4 ;
 *    	    X = 5 ;
 *    	    YY = 66661 ;
 *    	    XX = 66661 ;
 *    variables:
 *    	    int var(Y, X) ;
 *    	    float var_last(YY, XX) ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */

/* This program can also be used to test NetCDF.
 * Add #define TEST_NETCDF and compile with command:
 * gcc -I/netcdf/path/include last_large_var.c -o last_large_var -L/netcdf/path/lib -lnetcdf
 */
#ifdef TEST_NETCDF
#include <netcdf.h>
#include <netcdf_meta.h>
#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,nc_strerror(err)); \
    } \
}
#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %d but got %d\n", \
        __LINE__,__FILE__,exp, err); \
    } \
}
#define FileCreate(a,b,c,d,e)	nc_create(b,c,e)
#define DefDim			nc_def_dim
#define DefVar			nc_def_var
#define SetFill			nc_set_fill
#define ReDef			nc_redef
#define EndDef			nc_enddef
#define EndDef_			nc__enddef
#define FileClose		nc_close
#define MPI_Init(a,b)
#define MPI_Comm_rank(a,b)
#define MPI_Comm_size(a,b)
#define MPI_Finalize()
#define MPI_Bcast(a,b,c,d,e)
#else
#include <pnetcdf.h>
#include <testutils.h>
#define FileCreate	ncmpi_create
#define DefDim		ncmpi_def_dim
#define DefVar		ncmpi_def_var
#define SetFill		ncmpi_set_fill
#define ReDef		ncmpi_redef
#define EndDef		ncmpi_enddef
#define EndDef_		ncmpi__enddef
#define FileClose	ncmpi_close
#endif

static
int check_last_var(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new file ---------------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "X", 5, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "YY", 66661, &dimid[2]); CHECK_ERR
    err = DefDim(ncid, "XX", 66661, &dimid[3]); CHECK_ERR

    /* define only fixed-size variables and the last one is "big" */
    err = DefVar(ncid, "var", NC_INT, 1, dimid+1, &varid); CHECK_ERR
    err = DefVar(ncid, "var_last", NC_FLOAT, 2, dimid+2, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

    return nerrs;
}

static
int check_fix_var(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new CDF-1 file ----------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "X", 536870911, &dimid[0]); CHECK_ERR

    /* define only fixed-size variables and no one is "big"
     * make the starting offset of last one > 2GiB (illegal for CDF-1)
     */
    err = DefVar(ncid, "var1", NC_INT,   1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var2", NC_FLOAT, 1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var3", NC_SHORT, 1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var4", NC_INT,   1, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = FileClose(ncid); EXP_ERR(NC_EVARSIZE)

    /* create a new CDF-2 file ----------------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_OFFSET;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "X", 536870911, &dimid[0]); CHECK_ERR

    /* define only fixed-size variables and no one is "big"
     * make the starting offset of last one > 2GiB (legal for CDF-2)
     */
    err = DefVar(ncid, "var1", NC_INT,   1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var2", NC_FLOAT, 1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var3", NC_SHORT, 1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var4", NC_INT,   1, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

    return nerrs;
}

static
int check_fix_rec_var(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new file ---------------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "X", 5, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "YY", 66661, &dimid[2]); CHECK_ERR
    err = DefDim(ncid, "XX", 66661, &dimid[3]); CHECK_ERR

    /* define a record variable */
    err = DefVar(ncid, "var", NC_INT, 1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var_last", NC_FLOAT, 2, dimid+2, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); EXP_ERR(NC_EVARSIZE)
    err = FileClose(ncid); EXP_ERR(NC_EVARSIZE)

    return nerrs;
}

/* http://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#classic_format_limitations
 * If you use the unlimited dimension, record variables may exceed 2 GiB in
 * size, as long as the offset of the start of each record variable within a
 * record is less than 2 GiB - 4.
 */
static
int check_rec_var(const char *filename, int cmode, MPI_Info info)
{
    int err, nerrs=0, ncid, fill_mode, varid, dimid[3];

    /* create a new file ---------------------------------------------------*/
    cmode |= NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Z", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "Y", 66661,        &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "X", 66661,        &dimid[2]); CHECK_ERR

    /* define record variables: last one is large */
    err = DefVar(ncid, "var",       NC_INT,   1, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var_large", NC_FLOAT, 3, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

    /* create a new file ---------------------------------------------------*/
    cmode |= NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Z", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "Y", 1048576, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "X", 1000, &dimid[2]); CHECK_ERR

    /* define record variables: both starting offsets are < 2^31-4 */
    err = DefVar(ncid, "var1", NC_SHORT, 3, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var2", NC_SHORT, 3, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = FileClose(ncid); CHECK_ERR

    /* create a new file ---------------------------------------------------*/
    cmode |= NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Z", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "Y", 1048576, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "X", 1024, &dimid[2]); CHECK_ERR

    /* define record variables: some starting offsets are > 2^31-4 */
    err = DefVar(ncid, "var1", NC_SHORT, 3, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var2", NC_SHORT, 3, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var3", NC_SHORT, 3, dimid, &varid); CHECK_ERR
    err = DefVar(ncid, "var4", NC_SHORT, 3, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = FileClose(ncid);
    if (cmode & NC_64BIT_OFFSET || cmode & NC_64BIT_DATA) CHECK_ERR
    else EXP_ERR(NC_EVARSIZE)

    return nerrs;
}

/* http://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html#classic_format_limitations
 * If you don't use the unlimited dimension, only one variable can exceed 2 GiB
 * in size, but it can be as large as the underlying file system permits. It
 * must be the last variable in the dataset, and the offset to the beginning of
 * this variable must be less than about 2 GiB.
 */
static
int check_not_last_var(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new file ---------------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "X", 5, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "YY", 66661, &dimid[2]); CHECK_ERR
    err = DefDim(ncid, "XX", 66661, &dimid[3]); CHECK_ERR

    /* the large variable is not the last */
    err = DefVar(ncid, "var_large", NC_FLOAT, 2, dimid+2, &varid); CHECK_ERR
    err = DefVar(ncid, "var",       NC_INT,   1, dimid+1, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); EXP_ERR(NC_EVARSIZE)
    err = FileClose(ncid); EXP_ERR(NC_EVARSIZE)

    return nerrs;
}

static
int check_add_var(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new file ---------------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "X", 5, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "YY", 66661, &dimid[2]); CHECK_ERR
    err = DefDim(ncid, "XX", 66661, &dimid[3]); CHECK_ERR

    err = DefVar(ncid, "var", NC_INT, 1, dimid+1, &varid); CHECK_ERR
    err = DefVar(ncid, "var_last", NC_FLOAT, 2, dimid+2, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); CHECK_ERR

    /* add a new fixed-size variable */
    err = ReDef(ncid); CHECK_ERR
    err = DefVar(ncid, "var_new", NC_INT, 2, dimid, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    err = EndDef(ncid); EXP_ERR(NC_EVARSIZE)
    err = FileClose(ncid); EXP_ERR(NC_EVARSIZE)

    return nerrs;
}

static
int check_var_offset(const char *filename, MPI_Info info)
{
    int err, nerrs=0, ncid, cmode, fill_mode, varid, dimid[4];

    /* create a new file ---------------------------------------------------*/
    cmode = NC_CLOBBER;
    err = FileCreate(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = DefDim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = DefDim(ncid, "X", 5, &dimid[1]); CHECK_ERR
    err = DefDim(ncid, "YY", 66661, &dimid[2]); CHECK_ERR
    err = DefDim(ncid, "XX", 66661, &dimid[3]); CHECK_ERR

    err = DefVar(ncid, "var", NC_INT, 1, dimid+1, &varid); CHECK_ERR
    err = DefVar(ncid, "var_last", NC_FLOAT, 2, dimid+2, &varid); CHECK_ERR

    err = SetFill(ncid, NC_NOFILL, &fill_mode); CHECK_ERR
    /* make the file header size larger than 2 GiB */
    err = EndDef_(ncid, 2147483648LL, 1, 1, 1); EXP_ERR(NC_EVARSIZE)

    /* the above error keeps the program in define mode, thus close will
     * call enddef again and this time it will use default alignments, i.e.
     * ncmpi_enddef() is equivalent to ncmpi__enddef(ncid, 0, 1, 0, 1).
     * Thus, ncmpi_close() should return no error.
     */
    err = FileClose(ncid);
#if defined(TEST_NETCDF) && defined(NC_VERSION_MAJOR) && (NC_VERSION_MAJOR >= 4) && defined(NC_VERSION_MINOR) && (NC_VERSION_MINOR >= 7)
    EXP_ERR(NC_EVARSIZE) /* netCDF 4.6.0 and prior */
    /* See the following pull request in NetCDF
     * GitHub: https://github.com/Unidata/netcdf-c/pull/479
     */
#else
    CHECK_ERR
#endif

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,  /* ignored */
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int nerrs=0;

    nerrs += check_fix_var(out_path, info);
    nerrs += check_last_var(out_path, info);
    nerrs += check_fix_rec_var(out_path, info);

    nerrs += check_rec_var(out_path, 0, info);
    nerrs += check_rec_var(out_path, NC_64BIT_OFFSET, info);
    nerrs += check_rec_var(out_path, NC_64BIT_DATA, info);

    nerrs += check_not_last_var(out_path, info);
    nerrs += check_add_var(out_path, info);
    nerrs += check_var_offset(out_path, info);

#ifdef TEST_NETCDF
    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");
#endif

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {0};

    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "last large var in CDF-1/2", opt, test_io);

    MPI_Finalize();

    return err;
}
