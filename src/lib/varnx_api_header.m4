/* Begin of {put,get}_varn{kind} */

define(`VARNX',dnl
`dnl
int ncmpi_$1_varn$2(int ncid, int varid, int num, $3 void *buf,
              MPI_Offset bufcount, MPI_Datatype buftype);
')dnl
VARNX(get)
VARNX(get, _all)
VARNX(put,     , const)
VARNX(put, _all, const)

define(`VARN',dnl
`dnl
int ncmpi_$1_varn_$2$3(int ncid, int varid, int num, $4 *buf);
')dnl
VARN(get, text,         ,       char)
VARN(get, schar,        ,       signed char)
VARN(get, short,        ,       short)
VARN(get, int,          ,       int)
VARN(get, float,        ,       float)
VARN(get, double,       ,       double)
VARN(get, longlong,     ,       long long)
VARN(get, text,     _all,       char)
VARN(get, schar,    _all,       signed char)
VARN(get, short,    _all,       short)
VARN(get, int,      _all,       int)
VARN(get, float,    _all,       float)
VARN(get, double,   _all,       double)
VARN(get, longlong, _all,       long long)
VARN(put, text,         , const char)
VARN(put, schar,        , const signed char)
VARN(put, short,        , const short)
VARN(put, int,          , const int)
VARN(put, float,        , const float)
VARN(put, double,       , const double)
VARN(put, longlong,     , const long long)
VARN(put, text,     _all, const char)
VARN(put, schar,    _all, const signed char)
VARN(put, short,    _all, const short)
VARN(put, int,      _all, const int)
VARN(put, float,    _all, const float)
VARN(put, double,   _all, const double)
VARN(put, longlong, _all, const long long)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN(get, uchar,         ,       unsigned char)
VARN(get, ushort,        ,       unsigned short)
VARN(get, uint,          ,       unsigned int)
VARN(get, long,          ,       long)
VARN(get, ulonglong,     ,       unsigned long long)
VARN(get, uchar,     _all,       unsigned char)
VARN(get, ushort,    _all,       unsigned short)
VARN(get, uint,      _all,       unsigned int)
VARN(get, long,      _all,       long)
VARN(get, ulonglong, _all,       unsigned long long)
VARN(put, uchar,         , const unsigned char)
VARN(put, ushort,        , const unsigned short)
VARN(put, uint,          , const unsigned int)
VARN(put, long,          , const long)
VARN(put, ulonglong,     , const unsigned long long)
VARN(put, uchar,     _all, const unsigned char)
VARN(put, ushort,    _all, const unsigned short)
VARN(put, uint,      _all, const unsigned int)
VARN(put, long,      _all, const long)
VARN(put, ulonglong, _all, const unsigned long long)
/* End Skip Prototypes for Fortran binding */

define(`VARN1X',dnl
`dnl
int ncmpi_$1_varn1$2(int ncid, int varid, int num,
              MPI_Offset* const starts[], $3 void *buf,
              MPI_Offset bufcount, MPI_Datatype buftype);
')dnl
VARN1X(get)
VARN1X(get, _all)
VARN1X(put,     , const)
VARN1X(put, _all, const)

define(`VARN1',dnl
`dnl
int ncmpi_$1_varn1_$2$3(int ncid, int varid, int num,
              MPI_Offset* const starts[], $4 *buf);
')dnl
VARN1(get, text,         ,       char)
VARN1(get, schar,        ,       signed char)
VARN1(get, short,        ,       short)
VARN1(get, int,          ,       int)
VARN1(get, float,        ,       float)
VARN1(get, double,       ,       double)
VARN1(get, longlong,     ,       long long)
VARN1(get, text,     _all,       char)
VARN1(get, schar,    _all,       signed char)
VARN1(get, short,    _all,       short)
VARN1(get, int,      _all,       int)
VARN1(get, float,    _all,       float)
VARN1(get, double,   _all,       double)
VARN1(get, longlong, _all,       long long)
VARN1(put, text,         , const char)
VARN1(put, schar,        , const signed char)
VARN1(put, short,        , const short)
VARN1(put, int,          , const int)
VARN1(put, float,        , const float)
VARN1(put, double,       , const double)
VARN1(put, longlong,     , const long long)
VARN1(put, text,     _all, const char)
VARN1(put, schar,    _all, const signed char)
VARN1(put, short,    _all, const short)
VARN1(put, int,      _all, const int)
VARN1(put, float,    _all, const float)
VARN1(put, double,   _all, const double)
VARN1(put, longlong, _all, const long long)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARN1(get, uchar,         ,       unsigned char)
VARN1(get, ushort,        ,       unsigned short)
VARN1(get, uint,          ,       unsigned int)
VARN1(get, long,          ,       long)
VARN1(get, ulonglong,     ,       unsigned long long)
VARN1(get, uchar,     _all,       unsigned char)
VARN1(get, ushort,    _all,       unsigned short)
VARN1(get, uint,      _all,       unsigned int)
VARN1(get, long,      _all,       long)
VARN1(get, ulonglong, _all,       unsigned long long)
VARN1(put, uchar,         , const unsigned char)
VARN1(put, ushort,        , const unsigned short)
VARN1(put, uint,          , const unsigned int)
VARN1(put, long,          , const long)
VARN1(put, ulonglong,     , const unsigned long long)
VARN1(put, uchar,     _all, const unsigned char)
VARN1(put, ushort,    _all, const unsigned short)
VARN1(put, uint,      _all, const unsigned int)
VARN1(put, long,      _all, const long)
VARN1(put, ulonglong, _all, const unsigned long long)
/* End Skip Prototypes for Fortran binding */

define(`VARNAX',dnl
`dnl
int ncmpi_$1_varna$2(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              $3 void *buf,
              MPI_Offset bufcount, MPI_Datatype buftype);
')dnl
VARNAX(get)
VARNAX(get, _all)
VARNAX(put,     , const)
VARNAX(put, _all, const)

define(`VARNA',dnl
`dnl
int ncmpi_$1_varna_$2$3(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              $4 *buf);
')dnl
VARNA(get, text,         ,       char)
VARNA(get, schar,        ,       signed char)
VARNA(get, short,        ,       short)
VARNA(get, int,          ,       int)
VARNA(get, float,        ,       float)
VARNA(get, double,       ,       double)
VARNA(get, longlong,     ,       long long)
VARNA(get, text,     _all,       char)
VARNA(get, schar,    _all,       signed char)
VARNA(get, short,    _all,       short)
VARNA(get, int,      _all,       int)
VARNA(get, float,    _all,       float)
VARNA(get, double,   _all,       double)
VARNA(get, longlong, _all,       long long)
VARNA(put, text,         , const char)
VARNA(put, schar,        , const signed char)
VARNA(put, short,        , const short)
VARNA(put, int,          , const int)
VARNA(put, float,        , const float)
VARNA(put, double,       , const double)
VARNA(put, longlong,     , const long long)
VARNA(put, text,     _all, const char)
VARNA(put, schar,    _all, const signed char)
VARNA(put, short,    _all, const short)
VARNA(put, int,      _all, const int)
VARNA(put, float,    _all, const float)
VARNA(put, double,   _all, const double)
VARNA(put, longlong, _all, const long long)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARNA(get, uchar,         ,       unsigned char)
VARNA(get, ushort,        ,       unsigned short)
VARNA(get, uint,          ,       unsigned int)
VARNA(get, long,          ,       long)
VARNA(get, ulonglong,     ,       unsigned long long)
VARNA(get, uchar,     _all,       unsigned char)
VARNA(get, ushort,    _all,       unsigned short)
VARNA(get, uint,      _all,       unsigned int)
VARNA(get, long,      _all,       long)
VARNA(get, ulonglong, _all,       unsigned long long)
VARNA(put, uchar,         , const unsigned char)
VARNA(put, ushort,        , const unsigned short)
VARNA(put, uint,          , const unsigned int)
VARNA(put, long,          , const long)
VARNA(put, ulonglong,     , const unsigned long long)
VARNA(put, uchar,     _all, const unsigned char)
VARNA(put, ushort,    _all, const unsigned short)
VARNA(put, uint,      _all, const unsigned int)
VARNA(put, long,      _all, const long)
VARNA(put, ulonglong, _all, const unsigned long long)
/* End Skip Prototypes for Fortran binding */

define(`VARNSX',dnl
`dnl
int ncmpi_$1_varns$2(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              MPI_Offset* const strides[], $3 void *buf,
              MPI_Offset bufcount, MPI_Datatype buftype);
')dnl
VARNSX(get)
VARNSX(get, _all)
VARNSX(put,     , const)
VARNSX(put, _all, const)

define(`VARNS',dnl
`dnl
int ncmpi_$1_varns_$2$3(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              MPI_Offset* const strides[], $4 *buf);
')dnl
VARNS(get, text,         ,       char)
VARNS(get, schar,        ,       signed char)
VARNS(get, short,        ,       short)
VARNS(get, int,          ,       int)
VARNS(get, float,        ,       float)
VARNS(get, double,       ,       double)
VARNS(get, longlong,     ,       long long)
VARNS(get, text,     _all,       char)
VARNS(get, schar,    _all,       signed char)
VARNS(get, short,    _all,       short)
VARNS(get, int,      _all,       int)
VARNS(get, float,    _all,       float)
VARNS(get, double,   _all,       double)
VARNS(get, longlong, _all,       long long)
VARNS(put, text,         , const char)
VARNS(put, schar,        , const signed char)
VARNS(put, short,        , const short)
VARNS(put, int,          , const int)
VARNS(put, float,        , const float)
VARNS(put, double,       , const double)
VARNS(put, longlong,     , const long long)
VARNS(put, text,     _all, const char)
VARNS(put, schar,    _all, const signed char)
VARNS(put, short,    _all, const short)
VARNS(put, int,      _all, const int)
VARNS(put, float,    _all, const float)
VARNS(put, double,   _all, const double)
VARNS(put, longlong, _all, const long long)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARNS(get, uchar,         ,       unsigned char)
VARNS(get, ushort,        ,       unsigned short)
VARNS(get, uint,          ,       unsigned int)
VARNS(get, long,          ,       long)
VARNS(get, ulonglong,     ,       unsigned long long)
VARNS(get, uchar,     _all,       unsigned char)
VARNS(get, ushort,    _all,       unsigned short)
VARNS(get, uint,      _all,       unsigned int)
VARNS(get, long,      _all,       long)
VARNS(get, ulonglong, _all,       unsigned long long)
VARNS(put, uchar,         , const unsigned char)
VARNS(put, ushort,        , const unsigned short)
VARNS(put, uint,          , const unsigned int)
VARNS(put, long,          , const long)
VARNS(put, ulonglong,     , const unsigned long long)
VARNS(put, uchar,     _all, const unsigned char)
VARNS(put, ushort,    _all, const unsigned short)
VARNS(put, uint,      _all, const unsigned int)
VARNS(put, long,      _all, const long)
VARNS(put, ulonglong, _all, const unsigned long long)
/* End Skip Prototypes for Fortran binding */

define(`VARNMX',dnl
`dnl
int ncmpi_$1_varnm$2(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              MPI_Offset* const strides[], MPI_Offset* const imaps[],
              $3 void *buf,
              MPI_Offset bufcount, MPI_Datatype buftype);
')dnl
VARNMX(get)
VARNMX(get, _all)
VARNMX(put,     , const)
VARNMX(put, _all, const)

define(`VARNM',dnl
`dnl
int ncmpi_$1_varnm_$2$3(int ncid, int varid, int num,
              MPI_Offset* const starts[], MPI_Offset* const counts[],
              MPI_Offset* const strides[], MPI_Offset* const imaps[],
              $4 *buf);
')dnl
VARNM(get, text,         ,       char)
VARNM(get, schar,        ,       signed char)
VARNM(get, short,        ,       short)
VARNM(get, int,          ,       int)
VARNM(get, float,        ,       float)
VARNM(get, double,       ,       double)
VARNM(get, longlong,     ,       long long)
VARNM(get, text,     _all,       char)
VARNM(get, schar,    _all,       signed char)
VARNM(get, short,    _all,       short)
VARNM(get, int,      _all,       int)
VARNM(get, float,    _all,       float)
VARNM(get, double,   _all,       double)
VARNM(get, longlong, _all,       long long)
VARNM(put, text,         , const char)
VARNM(put, schar,        , const signed char)
VARNM(put, short,        , const short)
VARNM(put, int,          , const int)
VARNM(put, float,        , const float)
VARNM(put, double,       , const double)
VARNM(put, longlong,     , const long long)
VARNM(put, text,     _all, const char)
VARNM(put, schar,    _all, const signed char)
VARNM(put, short,    _all, const short)
VARNM(put, int,      _all, const int)
VARNM(put, float,    _all, const float)
VARNM(put, double,   _all, const double)
VARNM(put, longlong, _all, const long long)

/* Begin Skip Prototypes for Fortran binding */
/* skip types: uchar, ubyte, ushort, uint, long, ulonglong string */

VARNM(get, uchar,         ,       unsigned char)
VARNM(get, ushort,        ,       unsigned short)
VARNM(get, uint,          ,       unsigned int)
VARNM(get, long,          ,       long)
VARNM(get, ulonglong,     ,       unsigned long long)
VARNM(get, uchar,     _all,       unsigned char)
VARNM(get, ushort,    _all,       unsigned short)
VARNM(get, uint,      _all,       unsigned int)
VARNM(get, long,      _all,       long)
VARNM(get, ulonglong, _all,       unsigned long long)
VARNM(put, uchar,         , const unsigned char)
VARNM(put, ushort,        , const unsigned short)
VARNM(put, uint,          , const unsigned int)
VARNM(put, long,          , const long)
VARNM(put, ulonglong,     , const unsigned long long)
VARNM(put, uchar,     _all, const unsigned char)
VARNM(put, ushort,    _all, const unsigned short)
VARNM(put, uint,      _all, const unsigned int)
VARNM(put, long,      _all, const long)
VARNM(put, ulonglong, _all, const unsigned long long)
/* End Skip Prototypes for Fortran binding */

/* Begin of {put,get}_varn{kind} */

