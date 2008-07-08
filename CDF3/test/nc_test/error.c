/*********************************************************************
 *   Copyright 1996, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id: error.c 366 2004-05-26 21:43:25Z jianwli $
 *********************************************************************/

#include <stddef.h>	/* because gcc 2.7.2.2 doesn't define size_t */
			/* in <stdio.h> and it cannot hurt */
#include <stdio.h>
#include <stdarg.h>

#include <mpi.h>

extern int  nfails;             /* number of failures in specific test */
extern int  max_nmpt;		/* max. number of messages per test */

/* Prototypes */
void error(const char *fmt, ...);
void print(const char *fmt, ...);
int ifFail(const int expr, const int line, const char *file);
void print_n_size_t(size_t nelems, const MPI_Offset *array);

/*
 * Use for logging error conditions
 */
void
error(const char *fmt, ...)
{
    va_list args ;

    va_start(args, fmt) ;
    if (nfails <= max_nmpt)
	(void) vfprintf(stderr,fmt,args) ;
    va_end(args) ;
}

/*
 * Use for general printing (other than error conditions)
 * This also currently goes to stderr, but this could change
 */
void
print(const char *fmt, ...)
{
    va_list args ;

    va_start(args, fmt) ;
    (void) vfprintf(stderr,fmt,args) ;
    va_end(args) ;
}

/*
 * Called by macro IF
 */
int
ifFail(const int expr, const int line, const char *file)
{
    if (expr) {
	++nfails;
	error("\n\tFAILURE at line %d of %s: ", line, file);
    }
    return expr;
}

/* TODO:
 * This diagnostic doesn't fit very well with the diagnostic message 
 * "architecture" of this program.
 */
void
print_n_size_t(size_t nelems, const MPI_Offset *array)
{
	fprintf(stderr, "[");
	while(nelems-- > 0)
	{
		fprintf(stderr, "%ld", (long)*array++);
		if(nelems > 0)
			fprintf(stderr, " ");
	}
	fprintf(stderr, "]");
}
