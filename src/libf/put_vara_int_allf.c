/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface -infile=mpinetcdf.h -deffile=defs -debug
 * DO NOT EDIT
 */
#include "mpinetcdf_impl.h"


#ifdef F77_NAME_UPPER
#define ncmpif_put_vara_int_all_ NCMPIF_PUT_VARA_INT_ALL
#elif defined(F77_NAME_LOWER_2USCORE)
#define ncmpif_put_vara_int_all_ ncmpif_put_vara_int_all__
#elif !defined(F77_NAME_LOWER_USCORE)
#define ncmpif_put_vara_int_all_ ncmpif_put_vara_int_all
/* Else leave name alone */
#endif

FORTRAN_API void FORT_CALL ncmpif_put_vara_int_all_ ( int *v1, int *v2, int v3[], int v4[], int*v5, MPI_Fint *ierr )
{
    *ierr = ncmpi_put_vara_int_all( *v1, *v2, v3, v4, v5 );
}
