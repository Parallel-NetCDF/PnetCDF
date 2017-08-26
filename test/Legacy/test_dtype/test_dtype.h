#ifndef TEST_DTYPE_H
#define TEST_DTYPE_H

/* This is the only two macros you may need to edit for various etypes */
#define TEST_TYPE_INDEX 6 /* which type is used? [0..7] */
#define TEST_NCTYPE NC_DOUBLE
/* available memory types	VS	cooresponding MPI types 
	0 - unsigned char		MPI_UNSIGNED_CHAR
	1 - signed char			MPI_SIGNED_CHAR
	2 - char			MPI_CHAR
	3 - short			MPI_SHORT
	4 - int				MPI_INT
	5 - long			MPI_LONG
	6 - float			MPI_FLOAT
	7 - double 			MPI_DOUBLE
*/

#if TEST_TYPE_INDEX == 0
#define TEST_NATIVE_ETYPE unsigned char 
#define TEST_NATIVE_ETYPE_STR "unsigned char"
#elif TEST_TYPE_INDEX == 1
#define TEST_NATIVE_ETYPE signed char
#define TEST_NATIVE_ETYPE_STR "signed char"
#elif TEST_TYPE_INDEX == 2
#define TEST_NATIVE_ETYPE char
#define TEST_NATIVE_ETYPE_STR "char"
#elif TEST_TYPE_INDEX == 3
#define TEST_NATIVE_ETYPE short
#define TEST_NATIVE_ETYPE_STR "short"
#elif TEST_TYPE_INDEX == 4
#define TEST_NATIVE_ETYPE int
#define TEST_NATIVE_ETYPE_STR "int"
#elif TEST_TYPE_INDEX == 5
#define TEST_NATIVE_ETYPE long
#define TEST_NATIVE_ETYPE_STR "long"
#elif TEST_TYPE_INDEX == 6
#define TEST_NATIVE_ETYPE float
#define TEST_NATIVE_ETYPE_STR "float"
#elif TEST_TYPE_INDEX == 7
#define TEST_NATIVE_ETYPE double
#define TEST_NATIVE_ETYPE_STR "double"
#endif

#include <limits.h>
#ifdef INT_MAX
#define TEST_MAX_INT INT_MAX
#else
#define TEST_MAX_INT ( ~((int)(1) << (8*sizeof(int)-1)) )
#endif

#define TEST_SET_NCMPI_ETYPE(nc_etype, mpi_etype) 		\
{ 								\
  switch( (int)TEST_TYPE_INDEX ) {				\
  case 0:							\
    nc_etype = NC_BYTE; 					\
    mpi_etype = MPI_UNSIGNED_CHAR; 				\
    break;							\
  case 1:							\
    nc_etype = NC_BYTE; 					\
    mpi_etype = MPI_SIGNED_CHAR; 				\
    break;							\
  case 2:							\
    nc_etype = NC_CHAR; 					\
    mpi_etype = MPI_CHAR; 					\
    break;							\
  case 3:							\
    nc_etype = NC_SHORT; 					\
    mpi_etype = MPI_SHORT; 					\
    break;							\
  case 4:							\
    nc_etype = NC_INT; 						\
    mpi_etype = MPI_INT; 					\
    break;							\
  case 5:							\
    nc_etype = NC_INT; 						\
    mpi_etype = MPI_LONG; 					\
    break;							\
  case 6:							\
    nc_etype = NC_FLOAT; 					\
    mpi_etype = MPI_FLOAT; 					\
    break;							\
  case 7:							\
    nc_etype = NC_DOUBLE; 					\
    mpi_etype = MPI_DOUBLE; 					\
    break;							\
  }								\
}

#define TEST_GET_NCTYPE_STR(nc_etype)				\
(								\
  (nc_etype == NC_CHAR)?"NC_CHAR":(				\
  (nc_etype == NC_BYTE)?"NC_BYTE":(				\
  (nc_etype == NC_SHORT)?"NC_SHORT":(				\
  (nc_etype == NC_INT)?"NC_INT":(				\
  (nc_etype == NC_FLOAT)?"NC_FLOAT":(				\
  (nc_etype == NC_DOUBLE)?"NC_DOUBLE":"INVALID_TYPE")))))	\
)
    
#define TEST_NCTYPE_LEN(nc_etype)				\
(								\
  (nc_etype == NC_CHAR || nc_etype == NC_BYTE)	? 1 : (		\
  (nc_etype == NC_SHORT)			? 2 : (		\
  (nc_etype == NC_INT || nc_etype == NC_FLOAT)	? 4 : (		\
  (nc_etype == NC_DOUBLE) 			? 8 : 0 )))	\
)

#define TEST_HANDLE_ERR(status) 				\
{								\
  if ((status) != NC_NOERR) 					\
    printf( "%s\n", ncmpi_strerror((status)) );			\
}

#define TEST_REVERSE(array, nelems, type)			\
{								\
  int i;							\
  type tmp;							\
  for (i=0; i<(nelems)/2; i++) {				\
    tmp = (array)[i];						\
    (array)[i] = (array)[(nelems)-1-i];				\
    (array)[(nelems)-1-i] = tmp;				\
  }								\
}

#define TEST_PRINT_LIST(list, start, end, step)			\
{								\
  int i;							\
  for ( i=(start); (i-(start))*(i-(end))<=0; i+=(step) ) {	\
    if ( i==(start) )						\
      printf("%2d", (int)(list)[i]);				\
    else							\
      printf(", %2d", (int)(list)[i]);				\
  }								\
}

#define TEST_EXIT(err_code)					\
{								\
  MPI_Finalize();						\
  exit((err_code));						\
}

#endif
