/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2003 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef MPINETCDF_IMPL_H
#define MPINETCDF_IMPL_H
#include <stdlib.h>
#include "../lib/ncconfig.h"
#include "../lib/pnetcdf.h"

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
#define FORT_MIXED_LEN_DECL   , MPI_Fint
#define FORT_END_LEN_DECL
#define FORT_MIXED_LEN(a)     , MPI_Fint a
#define FORT_END_LEN(a)
#else
#define FORT_MIXED_LEN_DECL
#define FORT_END_LEN_DECL     , MPI_Fint
#define FORT_MIXED_LEN(a)
#define FORT_END_LEN(a)       , MPI_Fint a
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
int ncxVardim( int, int );

/* Define the internal values needed for Fortran support */

/* Fortran logicals */

/* Fortran logical values */

#endif




