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
#include <mpinetcdf_impl.h>


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
	return *value;
}

#ifdef F77_NAME_UPPER
#define nc_ignorefpe_ NC_IGNOREFPE
#elif defined(F77_NAME_LOWER_2USCORE)
#define nc_ignorefpe_  nc_ignorefpe__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nc_ignorefpe_  nc_ignorefpe
/* Else leave name alone */
#endif
#include <signal.h>
FORTRAN_API void FORT_CALL nc_ignorefpe_(int *doit)
{
	if(doit)
		(void) signal(SIGFPE, SIG_IGN);
}

#ifdef F77_NAME_UPPER
#define min_int_ MIN_INT
#elif defined(F77_NAME_LOWER_2USCORE)
#define min_int_  min_int__
#elif !defined(F77_NAME_LOWER_USCORE)
#define min_int_  min_int
/* Else leave name alone */
#endif
FORTRAN_API double min_int_() {
	return INT_MIN;
}

#ifdef F77_NAME_UPPER
#define max_int_ MAX_INT
#elif defined(F77_NAME_LOWER_2USCORE)
#define max_int_  max_int__
#elif !defined(F77_NAME_LOWER_USCORE)
#define max_int_  max_int
/* Else leave name alone */
#endif
FORTRAN_API double max_int_() {
	return INT_MAX;
}

#ifdef F77_NAME_UPPER
#define max_float_ MAX_FLOAT
#elif defined(F77_NAME_LOWER_2USCORE)
#define max_float_  max_float__
#elif !defined(F77_NAME_LOWER_USCORE)
#define max_float_  max_float
/* Else leave name alone */
#endif
FORTRAN_API double max_float_() {
	return FLT_MAX;
}

#ifdef F77_NAME_UPPER
#define max_double_ MAX_DOUBLE
#elif defined(F77_NAME_LOWER_2USCORE)
#define max_double_  max_double__
#elif !defined(F77_NAME_LOWER_USCORE)
#define max_double_  max_double
/* Else leave name alone */
#endif
FORTRAN_API double max_double_() {
	return LONG_MAX;
}
