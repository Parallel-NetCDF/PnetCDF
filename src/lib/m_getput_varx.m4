/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include "nc.h"
#include "ncx.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <assert.h>

#include "ncmpidtype.h"
#include "macro.h"


/* buffer layers:       
        
        User Level              buf     (user defined buffer of MPI_Datatype)
        MPI Datatype Level      cbuf    (contiguous buffer of ptype)
        NetCDF XDR Level        xbuf    (XDR I/O buffer)
*/

static int
ncmpii_mgetput_varm(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[], /* [num] */
                    int                rw_flag,     /* WRITE_REQ or READ_REQ */
                    int                io_method);  /* COLL_IO or INDEP_IO */

/*----< ncmpi_mput_var() >---------------------------------------------------*/
int
ncmpi_mput_var(int           ncid, 
               int           num, 
               int           varids[],    /* [num] */
               void         *bufs[],      /* [num] */
               MPI_Offset    bufcounts[], /* [num] */
               MPI_Datatype  datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_var_all() >-----------------------------------------------*/
int
ncmpi_mput_var_all(int           ncid, 
                   int           num, 
                   int           varids[],    /* [num] */
                   void         *bufs[],      /* [num] */
                   MPI_Offset    bufcounts[], /* [num] */
                   MPI_Datatype  datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

dnl
dnl MPUT_VAR(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MPUT_VAR',dnl
`dnl
/*----< ncmpi_mput_var_$1$5() >-----------------------------------------------*/
int
ncmpi_mput_var_$1$5(int ncid,
                    int num,
                    int varids[],
                    $2 *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              WRITE_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MPUT_VAR(text,      char,               MPI_CHAR,               INDEP_IO)
MPUT_VAR(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MPUT_VAR(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MPUT_VAR(short,     short,              MPI_SHORT,              INDEP_IO)
MPUT_VAR(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MPUT_VAR(int,       int,                MPI_INT,                INDEP_IO)
MPUT_VAR(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MPUT_VAR(long,      long,               MPI_LONG,               INDEP_IO)
MPUT_VAR(float,     float,              MPI_FLOAT,              INDEP_IO)
MPUT_VAR(double,    double,             MPI_DOUBLE,             INDEP_IO)
MPUT_VAR(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MPUT_VAR(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MPUT_VAR(text,      char,               MPI_CHAR,               COLL_IO, _all)
MPUT_VAR(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MPUT_VAR(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MPUT_VAR(short,     short,              MPI_SHORT,              COLL_IO, _all)
MPUT_VAR(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MPUT_VAR(int,       int,                MPI_INT,                COLL_IO, _all)
MPUT_VAR(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MPUT_VAR(long,      long,               MPI_LONG,               COLL_IO, _all)
MPUT_VAR(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MPUT_VAR(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MPUT_VAR(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MPUT_VAR(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mput_var1() >--------------------------------------------------*/
int
ncmpi_mput_var1(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_var1_all() >----------------------------------------------*/
int
ncmpi_mput_var1_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

dnl
dnl MPUT_VAR1(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MPUT_VAR1',dnl
`dnl
/*----< ncmpi_mput_var1_$1$5() >----------------------------------------------*/
int
ncmpi_mput_var1_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              WRITE_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MPUT_VAR1(text,      char,               MPI_CHAR,               INDEP_IO)
MPUT_VAR1(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MPUT_VAR1(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MPUT_VAR1(short,     short,              MPI_SHORT,              INDEP_IO)
MPUT_VAR1(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MPUT_VAR1(int,       int,                MPI_INT,                INDEP_IO)
MPUT_VAR1(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MPUT_VAR1(long,      long,               MPI_LONG,               INDEP_IO)
MPUT_VAR1(float,     float,              MPI_FLOAT,              INDEP_IO)
MPUT_VAR1(double,    double,             MPI_DOUBLE,             INDEP_IO)
MPUT_VAR1(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MPUT_VAR1(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MPUT_VAR1(text,      char,               MPI_CHAR,               COLL_IO, _all)
MPUT_VAR1(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MPUT_VAR1(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MPUT_VAR1(short,     short,              MPI_SHORT,              COLL_IO, _all)
MPUT_VAR1(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MPUT_VAR1(int,       int,                MPI_INT,                COLL_IO, _all)
MPUT_VAR1(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MPUT_VAR1(long,      long,               MPI_LONG,               COLL_IO, _all)
MPUT_VAR1(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MPUT_VAR1(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MPUT_VAR1(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MPUT_VAR1(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mput_vara() >--------------------------------------------------*/
int
ncmpi_mput_vara(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_vara_all() >----------------------------------------------*/
int
ncmpi_mput_vara_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

dnl
dnl MPUT_VARA(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MPUT_VARA',dnl
`dnl
/*----< ncmpi_mput_vara_$1$5() >----------------------------------------------*/
int
ncmpi_mput_vara_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              WRITE_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MPUT_VARA(text,      char,               MPI_CHAR,               INDEP_IO)
MPUT_VARA(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MPUT_VARA(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MPUT_VARA(short,     short,              MPI_SHORT,              INDEP_IO)
MPUT_VARA(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MPUT_VARA(int,       int,                MPI_INT,                INDEP_IO)
MPUT_VARA(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MPUT_VARA(long,      long,               MPI_LONG,               INDEP_IO)
MPUT_VARA(float,     float,              MPI_FLOAT,              INDEP_IO)
MPUT_VARA(double,    double,             MPI_DOUBLE,             INDEP_IO)
MPUT_VARA(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MPUT_VARA(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MPUT_VARA(text,      char,               MPI_CHAR,               COLL_IO, _all)
MPUT_VARA(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MPUT_VARA(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MPUT_VARA(short,     short,              MPI_SHORT,              COLL_IO, _all)
MPUT_VARA(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MPUT_VARA(int,       int,                MPI_INT,                COLL_IO, _all)
MPUT_VARA(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MPUT_VARA(long,      long,               MPI_LONG,               COLL_IO, _all)
MPUT_VARA(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MPUT_VARA(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MPUT_VARA(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MPUT_VARA(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mput_vars() >--------------------------------------------------*/
int
ncmpi_mput_vars(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                MPI_Offset* const  strides[],   /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_vars_all() >----------------------------------------------*/
int
ncmpi_mput_vars_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

dnl
dnl MPUT_VARS(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MPUT_VARS',dnl
`dnl
/*----< ncmpi_mput_vars_$1$5() >----------------------------------------------*/
int
ncmpi_mput_vars_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     MPI_Offset* const  strides[],   /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                              NULL, (void**)bufs, NULL, datatypes,
                              WRITE_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MPUT_VARS(text,      char,               MPI_CHAR,               INDEP_IO)
MPUT_VARS(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MPUT_VARS(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MPUT_VARS(short,     short,              MPI_SHORT,              INDEP_IO)
MPUT_VARS(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MPUT_VARS(int,       int,                MPI_INT,                INDEP_IO)
MPUT_VARS(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MPUT_VARS(long,      long,               MPI_LONG,               INDEP_IO)
MPUT_VARS(float,     float,              MPI_FLOAT,              INDEP_IO)
MPUT_VARS(double,    double,             MPI_DOUBLE,             INDEP_IO)
MPUT_VARS(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MPUT_VARS(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MPUT_VARS(text,      char,               MPI_CHAR,               COLL_IO, _all)
MPUT_VARS(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MPUT_VARS(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MPUT_VARS(short,     short,              MPI_SHORT,              COLL_IO, _all)
MPUT_VARS(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MPUT_VARS(int,       int,                MPI_INT,                COLL_IO, _all)
MPUT_VARS(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MPUT_VARS(long,      long,               MPI_LONG,               COLL_IO, _all)
MPUT_VARS(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MPUT_VARS(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MPUT_VARS(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MPUT_VARS(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mput_varm() >--------------------------------------------------*/
int
ncmpi_mput_varm(int                ncid, 
                int                num, 
                int                varids[],    /* [num] */
                MPI_Offset* const  starts[],    /* [num] */
                MPI_Offset* const  counts[],    /* [num] */
                MPI_Offset* const  strides[],   /* [num] */
                MPI_Offset* const  imaps[],     /* [num] */
                void              *bufs[],      /* [num] */
                MPI_Offset         bufcounts[], /* [num] */
                MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               WRITE_REQ, INDEP_IO);
}

/*----< ncmpi_mput_varm_all() >----------------------------------------------*/
int
ncmpi_mput_varm_all(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[]) /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               WRITE_REQ, COLL_IO);
}

dnl
dnl MPUT_VARM(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MPUT_VARM',dnl
`dnl
/*----< ncmpi_mput_varm_$1$5() >----------------------------------------------*/
int
ncmpi_mput_varm_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     MPI_Offset* const  strides[],   /* [num] */
                     MPI_Offset* const  imaps[],     /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                              imaps, (void**)bufs, NULL, datatypes,
                              WRITE_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MPUT_VARM(text,      char,               MPI_CHAR,               INDEP_IO)
MPUT_VARM(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MPUT_VARM(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MPUT_VARM(short,     short,              MPI_SHORT,              INDEP_IO)
MPUT_VARM(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MPUT_VARM(int,       int,                MPI_INT,                INDEP_IO)
MPUT_VARM(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MPUT_VARM(long,      long,               MPI_LONG,               INDEP_IO)
MPUT_VARM(float,     float,              MPI_FLOAT,              INDEP_IO)
MPUT_VARM(double,    double,             MPI_DOUBLE,             INDEP_IO)
MPUT_VARM(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MPUT_VARM(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MPUT_VARM(text,      char,               MPI_CHAR,               COLL_IO, _all)
MPUT_VARM(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MPUT_VARM(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MPUT_VARM(short,     short,              MPI_SHORT,              COLL_IO, _all)
MPUT_VARM(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MPUT_VARM(int,       int,                MPI_INT,                COLL_IO, _all)
MPUT_VARM(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MPUT_VARM(long,      long,               MPI_LONG,               COLL_IO, _all)
MPUT_VARM(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MPUT_VARM(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MPUT_VARM(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MPUT_VARM(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mget_var() >---------------------------------------------------*/
int
ncmpi_mget_var(int           ncid, 
               int           num, 
               int           varids[],     /* [num] */
               void         *bufs[],       /* [num] */
               MPI_Offset    bufcounts[],  /* [num] */
               MPI_Datatype  datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_var_all() >-----------------------------------------------*/
int
ncmpi_mget_var_all(int           ncid, 
                   int           num, 
                   int           varids[],     /* [num] */
                   void         *bufs[],       /* [num] */
                   MPI_Offset    bufcounts[],  /* [num] */
                   MPI_Datatype  datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

dnl
dnl MGET_VAR(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MGET_VAR',dnl
`dnl
/*----< ncmpi_mget_var_$1$5() >-----------------------------------------------*/
int
ncmpi_mget_var_$1$5(int ncid,
                    int num,
                    int varids[],
                    $2 *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, NULL, NULL, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              READ_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MGET_VAR(text,      char,               MPI_CHAR,               INDEP_IO)
MGET_VAR(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MGET_VAR(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MGET_VAR(short,     short,              MPI_SHORT,              INDEP_IO)
MGET_VAR(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MGET_VAR(int,       int,                MPI_INT,                INDEP_IO)
MGET_VAR(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MGET_VAR(long,      long,               MPI_LONG,               INDEP_IO)
MGET_VAR(float,     float,              MPI_FLOAT,              INDEP_IO)
MGET_VAR(double,    double,             MPI_DOUBLE,             INDEP_IO)
MGET_VAR(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MGET_VAR(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MGET_VAR(text,      char,               MPI_CHAR,               COLL_IO, _all)
MGET_VAR(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MGET_VAR(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MGET_VAR(short,     short,              MPI_SHORT,              COLL_IO, _all)
MGET_VAR(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MGET_VAR(int,       int,                MPI_INT,                COLL_IO, _all)
MGET_VAR(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MGET_VAR(long,      long,               MPI_LONG,               COLL_IO, _all)
MGET_VAR(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MGET_VAR(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MGET_VAR(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MGET_VAR(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mget_var1() >--------------------------------------------------*/
int
ncmpi_mget_var1(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_var1_all() >----------------------------------------------*/
int
ncmpi_mget_var1_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

dnl
dnl MGET_VAR1(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MGET_VAR1',dnl
`dnl
/*----< ncmpi_mget_var1_$1$5() >----------------------------------------------*/
int
ncmpi_mget_var1_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, NULL, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              READ_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MGET_VAR1(text,      char,               MPI_CHAR,               INDEP_IO)
MGET_VAR1(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MGET_VAR1(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MGET_VAR1(short,     short,              MPI_SHORT,              INDEP_IO)
MGET_VAR1(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MGET_VAR1(int,       int,                MPI_INT,                INDEP_IO)
MGET_VAR1(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MGET_VAR1(long,      long,               MPI_LONG,               INDEP_IO)
MGET_VAR1(float,     float,              MPI_FLOAT,              INDEP_IO)
MGET_VAR1(double,    double,             MPI_DOUBLE,             INDEP_IO)
MGET_VAR1(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MGET_VAR1(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MGET_VAR1(text,      char,               MPI_CHAR,               COLL_IO, _all)
MGET_VAR1(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MGET_VAR1(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MGET_VAR1(short,     short,              MPI_SHORT,              COLL_IO, _all)
MGET_VAR1(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MGET_VAR1(int,       int,                MPI_INT,                COLL_IO, _all)
MGET_VAR1(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MGET_VAR1(long,      long,               MPI_LONG,               COLL_IO, _all)
MGET_VAR1(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MGET_VAR1(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MGET_VAR1(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MGET_VAR1(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mget_vara() >--------------------------------------------------*/
int
ncmpi_mget_vara(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_vara_all() >----------------------------------------------*/
int
ncmpi_mget_vara_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

dnl
dnl MGET_VARA(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MGET_VARA',dnl
`dnl
/*----< ncmpi_mget_vara_$1$5() >----------------------------------------------*/
int
ncmpi_mget_vara_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, NULL,
                              NULL, (void**)bufs, NULL, datatypes,
                              READ_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MGET_VARA(text,      char,               MPI_CHAR,               INDEP_IO)
MGET_VARA(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MGET_VARA(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MGET_VARA(short,     short,              MPI_SHORT,              INDEP_IO)
MGET_VARA(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MGET_VARA(int,       int,                MPI_INT,                INDEP_IO)
MGET_VARA(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MGET_VARA(long,      long,               MPI_LONG,               INDEP_IO)
MGET_VARA(float,     float,              MPI_FLOAT,              INDEP_IO)
MGET_VARA(double,    double,             MPI_DOUBLE,             INDEP_IO)
MGET_VARA(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MGET_VARA(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MGET_VARA(text,      char,               MPI_CHAR,               COLL_IO, _all)
MGET_VARA(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MGET_VARA(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MGET_VARA(short,     short,              MPI_SHORT,              COLL_IO, _all)
MGET_VARA(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MGET_VARA(int,       int,                MPI_INT,                COLL_IO, _all)
MGET_VARA(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MGET_VARA(long,      long,               MPI_LONG,               COLL_IO, _all)
MGET_VARA(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MGET_VARA(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MGET_VARA(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MGET_VARA(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mget_vars() >--------------------------------------------------*/
int
ncmpi_mget_vars(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                MPI_Offset* const  strides[],    /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_vars_all() >----------------------------------------------*/
int
ncmpi_mget_vars_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    MPI_Offset* const  strides[],    /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               NULL, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

dnl
dnl MGET_VARS(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MGET_VARS',dnl
`dnl
/*----< ncmpi_mget_vars_$1$5() >----------------------------------------------*/
int
ncmpi_mget_vars_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     MPI_Offset* const  strides[],   /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                              NULL, (void**)bufs, NULL, datatypes,
                              READ_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MGET_VARS(text,      char,               MPI_CHAR,               INDEP_IO)
MGET_VARS(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MGET_VARS(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MGET_VARS(short,     short,              MPI_SHORT,              INDEP_IO)
MGET_VARS(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MGET_VARS(int,       int,                MPI_INT,                INDEP_IO)
MGET_VARS(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MGET_VARS(long,      long,               MPI_LONG,               INDEP_IO)
MGET_VARS(float,     float,              MPI_FLOAT,              INDEP_IO)
MGET_VARS(double,    double,             MPI_DOUBLE,             INDEP_IO)
MGET_VARS(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MGET_VARS(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MGET_VARS(text,      char,               MPI_CHAR,               COLL_IO, _all)
MGET_VARS(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MGET_VARS(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MGET_VARS(short,     short,              MPI_SHORT,              COLL_IO, _all)
MGET_VARS(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MGET_VARS(int,       int,                MPI_INT,                COLL_IO, _all)
MGET_VARS(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MGET_VARS(long,      long,               MPI_LONG,               COLL_IO, _all)
MGET_VARS(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MGET_VARS(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MGET_VARS(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MGET_VARS(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpi_mget_varm() >--------------------------------------------------*/
int
ncmpi_mget_varm(int                ncid, 
                int                num, 
                int                varids[],     /* [num] */
                MPI_Offset* const  starts[],     /* [num] */
                MPI_Offset* const  counts[],     /* [num] */
                MPI_Offset* const  strides[],    /* [num] */
                MPI_Offset* const  imaps[],      /* [num] */
                void              *bufs[],       /* [num] */
                MPI_Offset         bufcounts[],  /* [num] */
                MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               READ_REQ, INDEP_IO);
}

/*----< ncmpi_mget_varm_all() >----------------------------------------------*/
int
ncmpi_mget_varm_all(int                ncid, 
                    int                num, 
                    int                varids[],     /* [num] */
                    MPI_Offset* const  starts[],     /* [num] */
                    MPI_Offset* const  counts[],     /* [num] */
                    MPI_Offset* const  strides[],    /* [num] */
                    MPI_Offset* const  imaps[],      /* [num] */
                    void              *bufs[],       /* [num] */
                    MPI_Offset         bufcounts[],  /* [num] */
                    MPI_Datatype       datatypes[])  /* [num] */
{
    return ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                               imaps, bufs, bufcounts, datatypes,
                               READ_REQ, COLL_IO);
}

dnl
dnl MGET_VARM(btype_name, btype, mpi_type, mode, mode_name)
dnl
define(`MGET_VARM',dnl
`dnl
/*----< ncmpi_mget_varm_$1$5() >----------------------------------------------*/
int
ncmpi_mget_varm_$1$5(int                ncid,
                     int                num,
                     int                varids[],
                     MPI_Offset* const  starts[],    /* [num] */
                     MPI_Offset* const  counts[],    /* [num] */
                     MPI_Offset* const  strides[],   /* [num] */
                     MPI_Offset* const  imaps[],     /* [num] */
                     $2                *bufs[])
{
    int i, err;
    MPI_Datatype *datatypes;

    datatypes = (MPI_Datatype*) malloc(num * sizeof(MPI_Datatype));
    for (i=0; i<num; i++)
        datatypes[i] = $3;

    err = ncmpii_mgetput_varm(ncid, num, varids, starts, counts, strides,
                              imaps, (void**)bufs, NULL, datatypes,
                              READ_REQ, $4);
    free(datatypes);
    return err;
}
')dnl

MGET_VARM(text,      char,               MPI_CHAR,               INDEP_IO)
MGET_VARM(schar,     signed char,        MPI_BYTE,               INDEP_IO)
MGET_VARM(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      INDEP_IO)
MGET_VARM(short,     short,              MPI_SHORT,              INDEP_IO)
MGET_VARM(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     INDEP_IO)
MGET_VARM(int,       int,                MPI_INT,                INDEP_IO)
MGET_VARM(uint,      unsigned int,       MPI_UNSIGNED,           INDEP_IO)
MGET_VARM(long,      long,               MPI_LONG,               INDEP_IO)
MGET_VARM(float,     float,              MPI_FLOAT,              INDEP_IO)
MGET_VARM(double,    double,             MPI_DOUBLE,             INDEP_IO)
MGET_VARM(longlong,  long long,          MPI_LONG_LONG_INT,      INDEP_IO)
MGET_VARM(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, INDEP_IO)

MGET_VARM(text,      char,               MPI_CHAR,               COLL_IO, _all)
MGET_VARM(schar,     signed char,        MPI_BYTE,               COLL_IO, _all)
MGET_VARM(uchar,     unsigned char,      MPI_UNSIGNED_CHAR,      COLL_IO, _all)
MGET_VARM(short,     short,              MPI_SHORT,              COLL_IO, _all)
MGET_VARM(ushort,    unsigned short,     MPI_UNSIGNED_SHORT,     COLL_IO, _all)
MGET_VARM(int,       int,                MPI_INT,                COLL_IO, _all)
MGET_VARM(uint,      unsigned int,       MPI_UNSIGNED,           COLL_IO, _all)
MGET_VARM(long,      long,               MPI_LONG,               COLL_IO, _all)
MGET_VARM(float,     float,              MPI_FLOAT,              COLL_IO, _all)
MGET_VARM(double,    double,             MPI_DOUBLE,             COLL_IO, _all)
MGET_VARM(longlong,  long long,          MPI_LONG_LONG_INT,      COLL_IO, _all)
MGET_VARM(ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG, COLL_IO, _all)

/*----< ncmpii_mgetput_varm() >-----------------------------------------------*/
int
ncmpii_mgetput_varm(int                ncid, 
                    int                num, 
                    int                varids[],    /* [num] */
                    MPI_Offset* const  starts[],    /* [num] */
                    MPI_Offset* const  counts[],    /* [num] */
                    MPI_Offset* const  strides[],   /* [num] */
                    MPI_Offset* const  imaps[],     /* [num] */
                    void              *bufs[],      /* [num] */
                    MPI_Offset         bufcounts[], /* [num] */
                    MPI_Datatype       datatypes[], /* [num] */
                    int                rw_flag,     /* WRITE_REQ or READ_REQ */
                    int                io_method)   /* COLL_IO or INDEP_IO */
{
    int i, j, status=NC_NOERR, *req_ids=NULL, *statuses=NULL, min_st;
    NC *ncp=NULL;

    /* check if ncid is valid */
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        /* must return the error now, parallel program might hang */
        return status;

    /* check if it is in define mode */
    if (NC_indef(ncp)) status = NC_EINDEFINE;

    /* check file write permission if this is write request */
    if (status == NC_NOERR) {
        if (rw_flag == WRITE_REQ && NC_readonly(ncp)) status = NC_EPERM;
    }
    /* check whether collective or independent mode */
    if (status == NC_NOERR) {
        if (io_method == INDEP_IO)
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->independent_fh),
                                        MPI_COMM_SELF, 0);
        else if (io_method == COLL_IO)
            status = ncmpii_check_mpifh(ncp, &(ncp->nciop->collective_fh),
                                        ncp->nciop->comm, 1);
        /* else if (io_method == INDEP_COLL_IO) */
    }
    if (ncp->safe_mode == 1 && io_method == COLL_IO)
        MPI_Allreduce(&status, &min_st, 1, MPI_INT, MPI_MIN, ncp->nciop->comm);
    else
        min_st = status;

    if (min_st != NC_NOERR)
        return status;
  
    if (num > 0) {
        req_ids  = (int*) NCI_Malloc(2 * num * sizeof(int));
        statuses = req_ids + num;
    }

    /* for each request call ncmpi_igetput_varm() */
    for (i=0; i<num; i++) {
        NC_var *varp;
        MPI_Offset *start, *count, buflen;

        varp = ncmpii_NC_lookupvar(ncp, varids[i]);
        if (varp == NULL) continue; /* invalid varid, skip this request */

        /* when bufcounts == NULL, it means the same as counts[] */
        if (counts != NULL)
            for (buflen=1, j=0; j<varp->ndims; j++)
                 buflen *= counts[i][j];

        if (starts == NULL) {         /* var */
            GET_FULL_DIMENSIONS(start, count)
            for (buflen=1, j=0; j<varp->ndims; j++) buflen *= count[j];
            status = ncmpii_igetput_varm(ncp, varp, start, count, NULL,
                                         NULL, bufs[i], buflen,
                                         datatypes[i], &req_ids[i], rw_flag, 0);
            if (varp->ndims > 0) NCI_Free(start);
        } else if (counts == NULL) {  /* var1 */
            GET_FULL_DIMENSIONS(start, count)
            GET_ONE_COUNT(count)
            status = ncmpii_igetput_varm(ncp, varp, starts[i], count, NULL,
                                         NULL, bufs[i], 1,
                                         datatypes[i], &req_ids[i], rw_flag, 0);
            if (varp->ndims > 0) NCI_Free(count);
        } else if (strides == NULL) { /* vara */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i], NULL,
                                         NULL, bufs[i], buflen,
                                         datatypes[i], &req_ids[i], rw_flag, 0);
        } else if (imaps == NULL) {   /* vars */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], NULL, bufs[i], buflen,
                                         datatypes[i], &req_ids[i], rw_flag, 0);
        } else {                      /* varm */
            status = ncmpii_igetput_varm(ncp, varp, starts[i], counts[i],
                                         strides[i], imaps[i], bufs[i],
                                         buflen, datatypes[i],
                                         &req_ids[i], rw_flag, 0);
        }
    }

    if (status != NC_NOERR)
        return status;

    if (io_method == COLL_IO)
        status = ncmpi_wait_all(ncid, num, req_ids, statuses);
    else
        status = ncmpi_wait(ncid, num, req_ids, statuses);
    if (status != NC_NOERR)
        return status;

    if (num > 0)
        NCI_Free(req_ids);

    for (i=0; i<num; i++)
        if (statuses[i] != NC_NOERR)
            return statuses[i];

    return NC_NOERR;
}
