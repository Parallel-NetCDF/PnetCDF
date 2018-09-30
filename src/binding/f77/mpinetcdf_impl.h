/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *      See COPYRIGHT in top-level directory.
 */
/* $Id$ */

#ifndef _MPINETCDF_IMPL_H
#define _MPINETCDF_IMPL_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <string.h> /* memcpy(), memset() */
#include <pnetcdf.h>

/* Support Windows extension to specify calling convention */
#ifdef USE_FORT_STDCALL
#define FORT_CALL __stdcall
#elif defined (USE_FORT_CDECL)
#define FORT_CALL __cdecl
#else
#define FORT_CALL
#endif

/* Handle different mechanisms for passing Fortran CHARACTER to routines */
#ifdef USE_FORT_MIXED_STR_LEN
#define FORT_MIXED_LEN_DECL   , int
#define FORT_END_LEN_DECL
#define FORT_MIXED_LEN(a)     , int a
#define FORT_END_LEN(a)
#else
#define FORT_MIXED_LEN_DECL
#define FORT_END_LEN_DECL     , int
#define FORT_MIXED_LEN(a)
#define FORT_END_LEN(a)       , int a
#endif

/* Support Windows extension to specify which functions are exported from
   shared (DLL) libraries */
#ifdef HAVE_FORTRAN_API
# ifdef FORTRAN_EXPORTS
#  define FORTRAN_API __declspec(dllexport)
# else
#  define FORTRAN_API __declspec(dllimport)
# endif
#else
# define FORTRAN_API
#endif

/* Support an alternate approach for declaring a weak symbol supported by
   some versions of gcc */
#ifdef USE_WEAK_ATTRIBUTE
#define FUNC_ATTRIBUTES(name) __attribute__ ((weak,alias(#name)))
#else
#define FUNC_ATTRIBUTES(name)
#endif

/* Utility functions */
int ncmpixVardim( int, int );

extern FORTRAN_API int FORT_CALL nfmpi_xstrerror_ ( MPI_Fint *v1, char *v2 FORT_MIXED_LEN(d2) FORT_END_LEN(d2) );
extern FORTRAN_API int FORT_CALL nfmpi_xstrerrno_ ( MPI_Fint *v1, char *v2 FORT_MIXED_LEN(d2) FORT_END_LEN(d2) );
extern FORTRAN_API int FORT_CALL nfmpi_xinq_libvers_ ( char *v1 FORT_MIXED_LEN(d1) FORT_END_LEN(d1) );
extern FORTRAN_API int FORT_CALL nfmpi_issyserr_ ( MPI_Fint *v1 );
/* Define the internal values needed for Fortran support */

/* Fortran logicals */

/* Fortran logical values */

#endif

