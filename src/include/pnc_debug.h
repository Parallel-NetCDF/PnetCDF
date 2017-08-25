/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _PNC_DEBUG_H
#define _PNC_DEBUG_H

#include <stdio.h>

/* The two defines below are to be set manually.
#define PNETCDF_TRACE_MPI_COMM
#define PNETCDF_TRACE_MPI_IO
 */

/* C macros for TRACE MPI calls */
#ifdef PNETCDF_TRACE_MPI_COMM
#define TRACE_COMM(x) printf("TRACE-MPI-COMM: FILE %s FUNC %s() LINE %d calling %s()\n",__FILE__,__func__,__LINE__,#x),mpireturn=x
#else
#define TRACE_COMM(x) mpireturn=x
#endif

#ifdef PNETCDF_TRACE_MPI_IO
#define TRACE_IO(x) printf("TRACE-MPI-IO:   FILE %s FUNC %s() LINE %d calling %s()\n",__FILE__,__func__,__LINE__,#x),mpireturn=x
#else
#define TRACE_IO(x) mpireturn=x
#endif

#define CHECK_MPI_ERROR(mpi_errorcode, err_msg, nc_err) {                     \
    if (mpi_errorcode != MPI_SUCCESS) {                                       \
        char errorString[MPI_MAX_ERROR_STRING];                               \
        int rank, errorStringLen;                                             \
        MPI_Comm_rank(ncp->nciop->comm, &rank);                               \
        MPI_Error_string(mpi_errorcode, errorString, &errorStringLen);        \
        printf("%2d: MPI Failure at line %d of %s (%s : %s)\n",               \
               rank, __LINE__, __FILE__, err_msg, errorString);               \
        mpi_err = nc_err;                                                     \
    }                                                                         \
}

#ifdef PNETCDF_DEBUG
#define DEBUG_RETURN_ERROR(err) {                                       \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
    return err;                                                         \
}
#define DEBUG_ASSIGN_ERROR(status, err) {                               \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
    status = err;                                                       \
}
#define DEBUG_TRACE_ERROR(err) {                                        \
    char *_env_str = getenv("PNETCDF_VERBOSE_DEBUG_MODE");              \
    if (_env_str != NULL && *_env_str != '0') {                         \
        int _rank;                                                      \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                          \
        fprintf(stderr, "Rank %d: %s error at line %d of %s in %s\n",   \
        _rank,ncmpi_strerrno(err),__LINE__,__func__,__FILE__);          \
    }                                                                   \
}
#else
#define DEBUG_RETURN_ERROR(err) return err;
#define DEBUG_ASSIGN_ERROR(status, err) status = err;
#define DEBUG_TRACE_ERROR(err)
#endif

#endif /* _PNC_DEBUG_H */
