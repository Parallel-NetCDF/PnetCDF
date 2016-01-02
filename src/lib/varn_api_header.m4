/* Begin of {put,get}_varn{kind} */

define(`VARN_FLEXIBLE',dnl
`dnl
int ncmpi_$1_varn$4(int ncid, int varid, int num, MPI_Offset* const starts[],
              MPI_Offset* const counts[], $2 void *buf, MPI_Offset bufcount,
              MPI_Datatype buftype);
')dnl

dnl PnetCDF flexible APIs
VARN_FLEXIBLE(put, const, WRITE_REQ,     , INDEP_IO)
VARN_FLEXIBLE(put, const, WRITE_REQ, _all,  COLL_IO)
VARN_FLEXIBLE(get,      ,  READ_REQ,     , INDEP_IO)
VARN_FLEXIBLE(get,      ,  READ_REQ, _all,  COLL_IO)

define(`VARN',dnl
`dnl
int ncmpi_$1_varn_$6$4(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              $2 $7 *buf);
')dnl

VARN(put, const, WRITE_REQ,     , INDEP_IO, text,      char,               MPI_CHAR)
VARN(put, const, WRITE_REQ,     , INDEP_IO, schar,     signed char,        MPI_SIGNED_CHAR)
VARN(put, const, WRITE_REQ,     , INDEP_IO, short,     short,              MPI_SHORT)
VARN(put, const, WRITE_REQ,     , INDEP_IO, int,       int,                MPI_INT)
VARN(put, const, WRITE_REQ,     , INDEP_IO, float,     float,              MPI_FLOAT)
VARN(put, const, WRITE_REQ,     , INDEP_IO, double,    double,             MPI_DOUBLE)
VARN(put, const, WRITE_REQ,     , INDEP_IO, longlong,  long long,          MPI_LONG_LONG_INT)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN(put, const, WRITE_REQ,     , INDEP_IO, uchar,     unsigned char,      MPI_UNSIGNED_CHAR)
VARN(put, const, WRITE_REQ,     , INDEP_IO, ushort,    unsigned short,     MPI_UNSIGNED_SHORT)
VARN(put, const, WRITE_REQ,     , INDEP_IO, uint,      unsigned int,       MPI_UNSIGNED)
VARN(put, const, WRITE_REQ,     , INDEP_IO, long,      long,               MPI_LONG)
VARN(put, const, WRITE_REQ,     , INDEP_IO, ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
/* End Skip Prototypes for Fortran binding */

VARN(put, const, WRITE_REQ, _all,  COLL_IO, text,      char,               MPI_CHAR)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, schar,     signed char,        MPI_SIGNED_CHAR)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, short,     short,              MPI_SHORT)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, int,       int,                MPI_INT)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, float,     float,              MPI_FLOAT)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, double,    double,             MPI_DOUBLE)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, longlong,  long long,          MPI_LONG_LONG_INT)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN(put, const, WRITE_REQ, _all,  COLL_IO, uchar,     unsigned char,      MPI_UNSIGNED_CHAR)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, ushort,    unsigned short,     MPI_UNSIGNED_SHORT)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, uint,      unsigned int,       MPI_UNSIGNED)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, long,      long,               MPI_LONG)
VARN(put, const, WRITE_REQ, _all,  COLL_IO, ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
/* End Skip Prototypes for Fortran binding */

VARN(get,      ,  READ_REQ,     , INDEP_IO, text,      char,               MPI_CHAR)
VARN(get,      ,  READ_REQ,     , INDEP_IO, schar,     signed char,        MPI_SIGNED_CHAR)
VARN(get,      ,  READ_REQ,     , INDEP_IO, short,     short,              MPI_SHORT)
VARN(get,      ,  READ_REQ,     , INDEP_IO, int,       int,                MPI_INT)
VARN(get,      ,  READ_REQ,     , INDEP_IO, float,     float,              MPI_FLOAT)
VARN(get,      ,  READ_REQ,     , INDEP_IO, double,    double,             MPI_DOUBLE)
VARN(get,      ,  READ_REQ,     , INDEP_IO, longlong,  long long,          MPI_LONG_LONG_INT)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN(get,      ,  READ_REQ,     , INDEP_IO, uchar,     unsigned char,      MPI_UNSIGNED_CHAR)
VARN(get,      ,  READ_REQ,     , INDEP_IO, ushort,    unsigned short,     MPI_UNSIGNED_SHORT)
VARN(get,      ,  READ_REQ,     , INDEP_IO, uint,      unsigned int,       MPI_UNSIGNED)
VARN(get,      ,  READ_REQ,     , INDEP_IO, long,      long,               MPI_LONG)
VARN(get,      ,  READ_REQ,     , INDEP_IO, ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
/* End Skip Prototypes for Fortran binding */

VARN(get,      ,  READ_REQ, _all,  COLL_IO, text,      char,               MPI_CHAR)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, schar,     signed char,        MPI_SIGNED_CHAR)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, short,     short,              MPI_SHORT)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, int,       int,                MPI_INT)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, float,     float,              MPI_FLOAT)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, double,    double,             MPI_DOUBLE)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, longlong,  long long,          MPI_LONG_LONG_INT)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN(get,      ,  READ_REQ, _all,  COLL_IO, uchar,     unsigned char,      MPI_UNSIGNED_CHAR)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, ushort,    unsigned short,     MPI_UNSIGNED_SHORT)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, uint,      unsigned int,       MPI_UNSIGNED)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, long,      long,               MPI_LONG)
VARN(get,      ,  READ_REQ, _all,  COLL_IO, ulonglong, unsigned long long, MPI_UNSIGNED_LONG_LONG)
/* End Skip Prototypes for Fortran binding */

/* End of {put,get}_varn{kind} */

