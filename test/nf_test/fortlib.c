/*
 * $Id$
 *
 * This file contains support functions for FORTRAN code.  For example,
 * under HP-UX A.09.05, the U77 library doesn't contain the exit()
 * routine -- so we create one here.
 *
 * This file is heavily modified from nf_test/fortlib.c in the serial netcdf
 * source.
 */


#include <stdlib.h>
#include <limits.h>
#include <float.h>


#ifdef F77_NAME_UPPER
#define ud_exit_ UD_EXIT
#elif defined(F77_NAME_LOWER_2USCORE)
#define ud_exit_  ud_exit__
#elif !defined(F77_NAME_LOWER_USCORE)
#define ud_exit_  ud_exit
/* Else leave name alone */
#endif

FORTRAN_API void FORT_CALL ud_exit_(int *v1) {
	exit (*v1);
	return;
}

#ifdef F77_NAME_UPPER
#define ud_abort_ UD_ABORT
#elif defined(F77_NAME_LOWER_2USCORE)
#define ud_abort_  ud_abort__
#elif !defined(F77_NAME_LOWER_USCORE)
#define ud_abort_  ud_abort
/* Else leave name alone */
#endif

FORTRAN_API void FORT_CALL ud_abort_(int *v1) {
	abort();
	return;
}

#ifdef F77_NAME_UPPER
#define ud_rand_ UD_RAND
#elif defined(F77_NAME_LOWER_2USCORE)
#define ud_rand_  ud_rand__
#elif !defined(F77_NAME_LOWER_USCORE)
#define ud_rand_  ud_rand
/* Else leave name alone */
#endif

/*
* Return a pseudo-random value between 0.0 and 1.0.
*
* We don't use RAND_MAX here because not all compilation
* environments define it (e.g. gcc(1) under SunOS 4.1.3).
*/
FORTRAN_API double FORT_CALL ud_rand_(int *seed) {
	if (*seed != 0)
		srand(*seed);
	return (rand() % 32768 )/ 32767.0;
}

#ifdef F77_NAME_UPPER
#define ud_shift_ UD_SHIFT
#elif defined(F77_NAME_LOWER_2USCORE)
#define ud_shift_  ud_shift__
#elif !defined(F77_NAME_LOWER_USCORE)
#define ud_shift_  ud_shift
/* Else leave name alone */
#endif

FORTRAN_API int FORT_CALL ud_shift_(int * value, int *amount) {
	if(*amount < 0)
		*value >>= -*amount;
	else if (*amount > 0)
		*value <<= *amount;
	return value;
}
#ifdef F77_NAME_UPPER
#define nc_ignorefpe_ UD_SHIFT
#elif defined(F77_NAME_LOWER_2USCORE)
#define nc_ignorefpe_  ud_shift__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nc_ignorefpe_  ud_shift
/* Else leave name alone */
#endif
#include <signal.h>
FORTRAN_API void FORT_CALL nc_ignorefpe_(int *doit)
{
	if(doit)
		(void) signal(SIGFPE, SIG_IGN);
}
