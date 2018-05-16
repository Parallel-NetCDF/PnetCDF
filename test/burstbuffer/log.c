/*********************************************************************
 *
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#pragma warning(disable : 4996)

#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <stdint.h>
#include <string.h>
#include <pnetcdf.h>
#include <unistd.h>
#include <mpi.h>
#include <testutils.h>

#define MAXPROCESSES 1024
#define WIDTH 6
#define SINGLEPROCRANK 2
#define SINGLEPROCNP 5

int filecomp(char *a, char* b) {
    int ret;
    int cnt1 = 0;
    int cnt2 = 0;
    int flg = 0;

    FILE *fp1 ;
    FILE *fp2 ;

    fp1 = fopen(a, "r");
    if ( fp1 == NULL ) {
        printf("\n%s File can not be opened : \n", a);
        return -1;
    }

    fseek(fp1,0,SEEK_END);
    cnt1 = ftell(fp1);

    fp2 = fopen(b,"r");
    if ( fp2 == NULL ) {
        printf("\n%s File can not be opened : \n", b);
        return -1;
    }

    fseek(fp2,0,SEEK_END);
    cnt2 = ftell(fp2);

    fseek(fp1,0,SEEK_SET);
    fseek(fp2,0,SEEK_SET);

    if ( cnt1 != cnt2 ) {
        printf("\nFile contents are not same\n");
    }
    else{
        while( ! feof(fp1) ) {
            if ( fgetc(fp1) != fgetc(fp2) ) {
                flg = 1;
                break;
            }
        }

        if ( flg ) ret = 1;
        else ret = 0;
    }

    fclose(fp1);
    fclose(fp2);

    return ret;
}

/*
The test write a NP * NP matrix M, NP is the number of process:
put_vara:
Process N write N copy of it's rank to row N ([N, 0...WIDTH]) using different APIs on different variable
final result should be:
0 0 0 0 ...
1 1 1 1 ...
2 2 2 2 ...
.
.
.
*/
int simpletest(char* fname, int enable_log) {
    int buffer[MAXPROCESSES];
    MPI_Offset start[2], count[2];

    int i, j, ret, errlen;
    int NProc, MyRank, NP;      // Total process; Rank
    int fid;        // Data set ID
    int did[2];     // IDs of dimension
    int vid;        // IDs for variables
    int dims[2];
    char tmp[1024], tmp2[1024];
    MPI_Info Info;
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    if (NP == 1) {    // Act if there is WIDTH processes for easy debugging. Most debugger supports only single processes.
        NProc = SINGLEPROCNP;
        MyRank = SINGLEPROCRANK;
    }
    else{
        NProc = NP;
    }
    if (MyRank < MAXPROCESSES) {
        // Ensure each process have a independent buffer directory

        MPI_Info_create(&Info);
        if (enable_log) {
            MPI_Info_set(Info, "pnetcdf_log", "enable");
        }
        // Create new cdf file
        ret = ncmpi_create(MPI_COMM_WORLD, fname, NC_CLOBBER, Info, &fid);
        if (ret != NC_NOERR) {
            printf("Error create file\n");
            goto ERROR;
        }
        ret = ncmpi_set_fill(fid, NC_FILL, NULL);
        if (ret != NC_NOERR) {
            printf("Error set fill\n");
            goto ERROR;
        }
        ret = ncmpi_def_dim(fid, "X", NProc, did);  // X
        if (ret != NC_NOERR) {
            printf("Error def dim X\n");
            goto ERROR;
        }
        ret = ncmpi_def_dim(fid, "Y", NProc, did + 1);    // Y
        if (ret != NC_NOERR) {
            printf("Error def dim Y\n");
            goto ERROR;
        }
        ret = ncmpi_def_var(fid, "M", NC_INT, 2, did, vid);
        if (ret != NC_NOERR) {
            printf("Error def var M\n");
            goto ERROR;
        }
        ret = ncmpi_enddef(fid);
        if (ret != NC_NOERR) {
            printf("Error enddef\n");
            goto ERROR;
        }
        // Indep mode
        ret = ncmpi_begin_indep_data(fid);
        if (ret != NC_NOERR) {
            printf("Error begin indep\n");
            goto ERROR;
        }
        // We all write rank from now on
        for (i = 0; i < NProc; i++) {
            buffer[i] = MyRank;
        }

        // put_vara
        count[0] = 1;
        count[1] = NProc;
        start[0] = MyRank;
        start[1] = 0;
        ret = ncmpi_put_vara_int(fid, vid, start, count, buffer);
        if (ret != NC_NOERR) {
            MPI_Error_string(ret, tmp, &errlen);
            printf("Error put_varn: %d\n%s\n", errlen, tmp);
            goto ERROR;
        }
        // Collective mode
        ncmpi_end_indep_data(fid);
        if (ret != NC_NOERR) {
            printf("Error end indep");
            goto ERROR;
        }
        ncmpi_close(fid);       // Close file
        if (ret != NC_NOERR) {
            printf("Error close");
            goto ERROR;
        }
    }

ERROR:
    return 0;
}


/*
The demo write a NP * (NP * 4) matrix M, which is divided into 4 NP * WIDTH submatrixes, NP is the number of process:
There are 4 submatrix for demonstrating put_var1, put_vara, put_vars, put_varn, put_var is used in each variable to set it to 0 in the beginning
The same pattern is writen 4 time to 4 variables with the position of submatrix swapped:
0: put_var1, put_vara, put_vars, put_varn
1: put_vara, put_vars, put_varn, put_var1
2: put_vars, put_varn, put_var1, put_vara
3: put_varn, put_var1, put_vara, put_vars

put_var1:
Process N writes it's rank to [N, 0] leaving other cell unattained
final result should be:
0 0 0 0 ...
1 0 0 0 ...
2 0 0 0 ...
.
.
.
put_vara:
Process N write N copy of it's rank to row N ([N, 0...WIDTH]) using different APIs on different variable
final result should be:
0 0 0 0 ...
1 1 1 1 ...
2 2 2 2 ...
.
.
.
put_varas:
Process N write N copy of it's rank to every 2 columns on row N ([N, 0 2 4 ... WIDTH]) using different APIs on different variable
final result should be:
0 0 0 0 ...
1 0 1 0 ...
2 0 2 0 ...
.
.
.
put_varan:
Process N write N copy of it's rank to diagonal starts on row N ([(N ... N + WIDTH) % NP, 0 ... WIDTH]) using different APIs on different variable
final result should be:
0 0 0 0 ...
1 0 0 0 ...
2 1 0 0 ...
3 2 1 0 ...
.
.
*/

int test(char* fname, int enable_log) {
    int buffer[MAXPROCESSES];
    MPI_Offset start[MAXPROCESSES][2], count[MAXPROCESSES][2];
    MPI_Offset *sp[MAXPROCESSES], *cp[MAXPROCESSES];
    MPI_Offset stride[2];

    int i, j, ret;
    int NProc, MyRank, NP;      // Total process; Rank
    int fid;        // Data set ID
    int did[2];     // IDs of dimension
    int vid[4];        // IDs for variables
    int dims[2];
    char tmp[1024];
    MPI_Info Info;

    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    if (NP == 1) {    // Act if there is WIDTH processes for easy debugging. Most debugger supports only single proccesses.
        NProc = SINGLEPROCNP;
        MyRank = SINGLEPROCRANK;
    }
    else{
        NProc = NP;
    }

    if (MyRank < MAXPROCESSES) {
        // Ensure each process have a independent buffer directory

        MPI_Info_create(&Info);
        if (enable_log) {
            MPI_Info_set(Info, "pnetcdf_log", "enable");
        }

        // Create new cdf file
        ret = ncmpi_create(MPI_COMM_WORLD, fname, NC_CLOBBER, Info, &fid);
        if (ret != NC_NOERR) {
            printf("Error create file\n");
            goto ERROR;
        }
        ret = ncmpi_set_fill(fid, NC_FILL, NULL);
        if (ret != NC_NOERR) {
            printf("Error set fill\n");
            goto ERROR;
        }
        ret = ncmpi_def_dim(fid, "X", NProc, did);  // X
        if (ret != NC_NOERR) {
            printf("Error def dim X\n");
            goto ERROR;
        }
        ret = ncmpi_def_dim(fid, "Y", NProc * 4, did + 1);    // Y
        if (ret != NC_NOERR) {
            printf("Error def dim Y\n");
            goto ERROR;
        }
        ret = ncmpi_def_var(fid, "M0", NC_INT, 2, did, vid + 0);
        if (ret != NC_NOERR) {
            printf("Error def var M0\n");
            goto ERROR;
        }
        ret = ncmpi_def_var(fid, "M1", NC_INT, 2, did, vid + 1);
        if (ret != NC_NOERR) {
            printf("Error def var M1\n");
            goto ERROR;
        }
        ret = ncmpi_def_var(fid, "M2", NC_INT, 2, did, vid + 2);
        if (ret != NC_NOERR) {
            printf("Error def var M2\n");
            goto ERROR;
        }
        ret = ncmpi_def_var(fid, "M3", NC_INT, 2, did, vid + 3);
        if (ret != NC_NOERR) {
            printf("Error def var M3\n");
            goto ERROR;
        }
        ret = ncmpi_enddef(fid);
        if (ret != NC_NOERR) {
            printf("Error enddef\n");
            goto ERROR;
        }

        // We all write rank from now on
        for (i = 0; i < NProc; i++) {
            buffer[i] = MyRank;
        }

        // put_var1
        for (i = 0; i < 4; i++) {
            for (j = 0; j < NProc; j++) {
                start[0][0] = MyRank;
                start[0][1] = i * NProc + j;
                ret = ncmpi_put_var1_int_all(fid, vid[i], start[0], buffer);
                if (ret != NC_NOERR) {
                    printf("Error put_var1\n");
                    goto ERROR;
                }
            }
        }

        // put_vara
        for (i = 0; i < 4; i++) {
            start[0][0] = 0;
            start[0][1] = ((i + 1) % 4) * NProc + MyRank;
            count[0][0] = NProc;
            count[0][1] = 1;
            ret = ncmpi_put_vara_int_all(fid, vid[i], start[0], count[0], buffer);
            if (ret != NC_NOERR) {
                printf("Error put_vara\n");
                goto ERROR;
            }
        }

        // put_vars
        for (i = 0; i < 4; i++) {
            start[0][0] = MyRank;
            start[0][1] = ((i + 2) % 4) * NProc + (MyRank % 2);
            count[0][0] = 1;
            count[0][1] = NProc / 2;
            stride[0] = 1;
            stride[1] = 2;
            ret = ncmpi_put_vars_int_all(fid, vid[i], start[0], count[0], stride, buffer);
            if (ret != NC_NOERR) {
                printf("Error put_vars\n");
                goto ERROR;
            }
        }

        // put_varn
        for (j = 0; j < 4; j++) {
            for (i = 0; i < NProc; i++) {
                count[i][0] = 1;
                count[i][1] = 1;
                start[i][0] = (MyRank + i) % NProc;
                start[i][1] = i + ((j + 3) % 4) * NProc;
                sp[i] = (MPI_Offset*)start[i];
                cp[i] = (MPI_Offset*)count[i];
            }
            ret = ncmpi_put_varn_int_all(fid, vid[j], NProc, sp, cp, buffer);
            if (ret != NC_NOERR) {
                printf("Error put_varn\n");
                goto ERROR;
            }
        }

        // Commit log into cdf file

        ret = ncmpi_close(fid);       // Close file
        if (ret != NC_NOERR) {
            printf("Error close");
            goto ERROR;
        }
    }

ERROR:;
    return 0;
}

int main(int argc, char* argv[]) {
    int err = 0;
    int rank;

    if (argc < 2) {
        printf("./log <output>\n");
        return 0;
    }

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for checking log functionality", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    simpletest("test.cdf", 0);
    simpletest("test_log.cdf", 1);

    if (filecomp("test.cdf", "test_log.cdf") != 0) {
        printf("Error in %s line %d: Flushed result mismatch\n", __FILE__,__LINE__);
        err += 1;
    }
    test("test.cdf", 0);
    test("test_log.cdf", 1);

    if (filecomp("test.cdf", "test_log.cdf") != 0) {
        printf("Error in %s line %d: Flushed result mismatch\n", __FILE__,__LINE__);
        err += 1;
    }

    if (err == 0) {
        printf("pass\n");
    }
    else{
        printf("fail with %d mismatches\n", err);
    }

    system("rm -f test.cdf");
    system("rm -f test_log.cdf");

ERROR:;
    MPI_Finalize();
    return err;
}
