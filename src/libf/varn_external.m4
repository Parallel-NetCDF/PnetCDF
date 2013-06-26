!
! Begin of varn subroutines
!

define(`VARNX',dnl
`dnl
      integer nfmpi_$1_varn$3
!     INTEGER FUNCTION nfmpi_$1_varn$3(ncid, varid, num, buf, bufcount, buftype)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      <type>,                        INTENT($2)  :: buf(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                      INTEGER,                       INTENT(IN)  :: buftype
!     END FUNCTION nfmpi_$1_varn$3
      external nfmpi_$1_varn$3
')dnl
!
! flexible APIs
!
VARNX(get, OUT)
VARNX(get, OUT, _all)
VARNX(put, IN)
VARNX(put, IN,  _all)

define(`VARN',dnl
`dnl
      integer nfmpi_$1_varn_$3$4
!     INTEGER FUNCTION nfmpi_$1_varn_$3$4(ncid, varid, num, $5)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      $6,			    INTENT($2)  :: $7
!     END FUNCTION nfmpi_$1_varn_$3$4
      external nfmpi_$1_varn_$3$4
')dnl
VARN(get, OUT, text,       , text,   CHARACTER(len=*), text)
VARN(get, OUT, int1,       , i1vals, INTEGER*1,        i1vals(*))
VARN(get, OUT, int2,       , i2vals, INTEGER*2,        i2vals(*))
VARN(get, OUT, int,        , ivals,  INTEGER,          ivals(*))
VARN(get, OUT, real,       , rvals,  REAL,             rvals(*))
VARN(get, OUT, double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARN(get, OUT, int8,       , i8vals, INTEGER*8,        i8vals(*))
VARN(get, OUT, text,   _all, text,   CHARACTER(len=*), text)
VARN(get, OUT, int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARN(get, OUT, int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARN(get, OUT, int,    _all, ivals,  INTEGER,          ivals(*))
VARN(get, OUT, real,   _all, rvals,  REAL,             rvals(*))
VARN(get, OUT, double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARN(get, OUT, int8,   _all, i8vals, INTEGER*8,        i8vals(*))

VARN(put, IN,  text,       , text,   CHARACTER(len=*), text)
VARN(put, IN,  int1,       , i1vals, INTEGER*1,        i1vals(*))
VARN(put, IN,  int2,       , i2vals, INTEGER*2,        i2vals(*))
VARN(put, IN,  int,        , ivals,  INTEGER,          ivals(*))
VARN(put, IN,  real,       , rvals,  REAL,             rvals(*))
VARN(put, IN,  double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARN(put, IN,  int8,       , i8vals, INTEGER*8,        i8vals(*))
VARN(put, IN,  text,   _all, text,   CHARACTER(len=*), text)
VARN(put, IN,  int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARN(put, IN,  int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARN(put, IN,  int,    _all, ivals,  INTEGER,          ivals(*))
VARN(put, IN,  real,   _all, rvals,  REAL,             rvals(*))
VARN(put, IN,  double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARN(put, IN,  int8,   _all, i8vals, INTEGER*8,        i8vals(*))

define(`VARN1X',dnl
`dnl
      integer nfmpi_$1_varn1$3
!     INTEGER FUNCTION nfmpi_$1_varn1$3(ncid, varid, num, starts, buf, bufcount, buftype)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      <type>,                        INTENT($2)  :: buf(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                      INTEGER,                       INTENT(IN)  :: buftype
!     END FUNCTION nfmpi_$1_varn1$3
      external nfmpi_$1_varn1$3
')dnl
!
! flexible APIs
!
VARN1X(get, OUT)
VARN1X(get, OUT, _all)
VARN1X(put, IN)
VARN1X(put, IN,  _all)

define(`VARN1',dnl
`dnl
      integer nfmpi_$1_varn1_$3$4
!     INTEGER FUNCTION nfmpi_$1_varn1_$3$4(ncid, varid, num, starts, $5)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      $6,			    INTENT($2)  :: $7
!     END FUNCTION nfmpi_$1_varn1_$3$4
      external nfmpi_$1_varn1_$3$4
')dnl
VARN1(get, OUT, text,       , text,   CHARACTER,        text)
VARN1(get, OUT, int1,       , i1vals, INTEGER*1,        i1vals(*))
VARN1(get, OUT, int2,       , i2vals, INTEGER*2,        i2vals(*))
VARN1(get, OUT, int,        , ivals,  INTEGER,          ivals(*))
VARN1(get, OUT, real,       , rvals,  REAL,             rvals(*))
VARN1(get, OUT, double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARN1(get, OUT, int8,       , i8vals, INTEGER*8,        i8vals(*))
VARN1(get, OUT, text,   _all, text,   CHARACTER,        text)
VARN1(get, OUT, int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARN1(get, OUT, int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARN1(get, OUT, int,    _all, ivals,  INTEGER,          ivals(*))
VARN1(get, OUT, real,   _all, rvals,  REAL,             rvals(*))
VARN1(get, OUT, double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARN1(get, OUT, int8,   _all, i8vals, INTEGER*8,        i8vals(*))

VARN1(put, IN,  text,       , text,   CHARACTER,        text)
VARN1(put, IN,  int1,       , i1vals, INTEGER*1,        i1vals(*))
VARN1(put, IN,  int2,       , i2vals, INTEGER*2,        i2vals(*))
VARN1(put, IN,  int,        , ivals,  INTEGER,          ivals(*))
VARN1(put, IN,  real,       , rvals,  REAL,             rvals(*))
VARN1(put, IN,  double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARN1(put, IN,  int8,       , i8vals, INTEGER*8,        i8vals(*))
VARN1(put, IN,  text,   _all, text,   CHARACTER,        text)
VARN1(put, IN,  int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARN1(put, IN,  int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARN1(put, IN,  int,    _all, ivals,  INTEGER,          ivals(*))
VARN1(put, IN,  real,   _all, rvals,  REAL,             rvals(*))
VARN1(put, IN,  double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARN1(put, IN,  int8,   _all, i8vals, INTEGER*8,        i8vals(*))

define(`VARNAX',dnl
`dnl
      integer nfmpi_$1_varna$3
!     INTEGER FUNCTION nfmpi_$1_varna$3(ncid, varid, num, starts, counts, buf, bufcount, buftype)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      <type>,                        INTENT($2)  :: buf(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                      INTEGER,                       INTENT(IN)  :: buftype
!     END FUNCTION nfmpi_$1_varna$3
      external nfmpi_$1_varna$3
')dnl
!
! flexible APIs
!
VARNAX(get, OUT)
VARNAX(get, OUT, _all)
VARNAX(put, IN)
VARNAX(put, IN,  _all)

define(`VARNA',dnl
`dnl
      integer nfmpi_$1_varna_$3$4
!     INTEGER FUNCTION nfmpi_$1_varna_$3$4(ncid, varid, num, starts, counts, $5)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      $6,			    INTENT($2)  :: $7
!     END FUNCTION nfmpi_$1_varna_$3$4
      external nfmpi_$1_varna_$3$4
')dnl
VARNA(get, OUT, text,       , text,   CHARACTER(len=*), text)
VARNA(get, OUT, int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNA(get, OUT, int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNA(get, OUT, int,        , ivals,  INTEGER,          ivals(*))
VARNA(get, OUT, real,       , rvals,  REAL,             rvals(*))
VARNA(get, OUT, double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNA(get, OUT, int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNA(get, OUT, text,   _all, text,   CHARACTER(len=*), text)
VARNA(get, OUT, int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNA(get, OUT, int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNA(get, OUT, int,    _all, ivals,  INTEGER,          ivals(*))
VARNA(get, OUT, real,   _all, rvals,  REAL,             rvals(*))
VARNA(get, OUT, double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNA(get, OUT, int8,   _all, i8vals, INTEGER*8,        i8vals(*))

VARNA(put, IN,  text,       , text,   CHARACTER(len=*), text)
VARNA(put, IN,  int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNA(put, IN,  int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNA(put, IN,  int,        , ivals,  INTEGER,          ivals(*))
VARNA(put, IN,  real,       , rvals,  REAL,             rvals(*))
VARNA(put, IN,  double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNA(put, IN,  int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNA(put, IN,  text,   _all, text,   CHARACTER(len=*), text)
VARNA(put, IN,  int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNA(put, IN,  int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNA(put, IN,  int,    _all, ivals,  INTEGER,          ivals(*))
VARNA(put, IN,  real,   _all, rvals,  REAL,             rvals(*))
VARNA(put, IN,  double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNA(put, IN,  int8,   _all, i8vals, INTEGER*8,        i8vals(*))

define(`VARNSX',dnl
`dnl
      integer nfmpi_$1_varns$3
!     INTEGER FUNCTION nfmpi_$1_varns$3(ncid, varid, num, starts, counts, strides, buf, bufcount, buftype)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: strides(*)
!                      <type>,                        INTENT($2)  :: buf(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                      INTEGER,                       INTENT(IN)  :: buftype
!     END FUNCTION nfmpi_$1_varns$3
      external nfmpi_$1_varns$3
')dnl
!
! flexible APIs
!
VARNSX(get, OUT)
VARNSX(get, OUT, _all)
VARNSX(put, IN)
VARNSX(put, IN,  _all)

define(`VARNS',dnl
`dnl
      integer nfmpi_$1_varns_$3$4
!     INTEGER FUNCTION nfmpi_$1_varns_$3$4(ncid, varid, num, starts, counts, strides, $5)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: strides(*)
!                      $6,			    INTENT($2)  :: $7
!     END FUNCTION nfmpi_$1_varns_$3$4
      external nfmpi_$1_varns_$3$4
')dnl
VARNS(get, OUT, text,       , text,   CHARACTER(len=*), text)
VARNS(get, OUT, int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNS(get, OUT, int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNS(get, OUT, int,        , ivals,  INTEGER,          ivals(*))
VARNS(get, OUT, real,       , rvals,  REAL,             rvals(*))
VARNS(get, OUT, double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNS(get, OUT, int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNS(get, OUT, text,   _all, text,   CHARACTER(len=*), text)
VARNS(get, OUT, int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNS(get, OUT, int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNS(get, OUT, int,    _all, ivals,  INTEGER,          ivals(*))
VARNS(get, OUT, real,   _all, rvals,  REAL,             rvals(*))
VARNS(get, OUT, double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNS(get, OUT, int8,   _all, i8vals, INTEGER*8,        i8vals(*))

VARNS(put, IN,  text,       , text,   CHARACTER(len=*), text)
VARNS(put, IN,  int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNS(put, IN,  int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNS(put, IN,  int,        , ivals,  INTEGER,          ivals(*))
VARNS(put, IN,  real,       , rvals,  REAL,             rvals(*))
VARNS(put, IN,  double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNS(put, IN,  int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNS(put, IN,  text,   _all, text,   CHARACTER(len=*), text)
VARNS(put, IN,  int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNS(put, IN,  int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNS(put, IN,  int,    _all, ivals,  INTEGER,          ivals(*))
VARNS(put, IN,  real,   _all, rvals,  REAL,             rvals(*))
VARNS(put, IN,  double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNS(put, IN,  int8,   _all, i8vals, INTEGER*8,        i8vals(*))

define(`VARNMX',dnl
`dnl
      integer nfmpi_$1_varnm$3
!     INTEGER FUNCTION nfmpi_$1_varnm$3(ncid, varid, num, starts, counts, strides, imaps, buf, bufcount, buftype)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: strides(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imaps(*)
!                      <type>,                        INTENT($2)  :: buf(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                      INTEGER,                       INTENT(IN)  :: buftype
!     END FUNCTION nfmpi_$1_varnm$3
      external nfmpi_$1_varnm$3
')dnl
!
! flexible APIs
!
VARNMX(get, OUT)
VARNMX(get, OUT, _all)
VARNMX(put, IN)
VARNMX(put, IN,  _all)

define(`VARNM',dnl
`dnl
      integer nfmpi_$1_varnm_$3$4
!     INTEGER FUNCTION nfmpi_$1_varnm_$3$4(ncid, varid, num, starts, counts, strides, imaps, $5)
!                      INTEGER,                       INTENT(IN)  :: ncid
!                      INTEGER,                       INTENT(IN)  :: varid
!                      INTEGER,                       INTENT(IN)  :: num
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: strides(*)
!                      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imaps(*)
!                      $6,			    INTENT($2)  :: $7
!     END FUNCTION nfmpi_$1_varnm_$3$4
      external nfmpi_$1_varnm_$3$4
')dnl
VARNM(get, OUT, text,       , text,   CHARACTER(len=*), text)
VARNM(get, OUT, int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNM(get, OUT, int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNM(get, OUT, int,        , ivals,  INTEGER,          ivals(*))
VARNM(get, OUT, real,       , rvals,  REAL,             rvals(*))
VARNM(get, OUT, double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNM(get, OUT, int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNM(get, OUT, text,   _all, text,   CHARACTER(len=*), text)
VARNM(get, OUT, int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNM(get, OUT, int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNM(get, OUT, int,    _all, ivals,  INTEGER,          ivals(*))
VARNM(get, OUT, real,   _all, rvals,  REAL,             rvals(*))
VARNM(get, OUT, double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNM(get, OUT, int8,   _all, i8vals, INTEGER*8,        i8vals(*))

VARNM(put, IN,  text,       , text,   CHARACTER(len=*), text)
VARNM(put, IN,  int1,       , i1vals, INTEGER*1,        i1vals(*))
VARNM(put, IN,  int2,       , i2vals, INTEGER*2,        i2vals(*))
VARNM(put, IN,  int,        , ivals,  INTEGER,          ivals(*))
VARNM(put, IN,  real,       , rvals,  REAL,             rvals(*))
VARNM(put, IN,  double,     , dvals,  DOUBLE PRECISION, dvals(*))
VARNM(put, IN,  int8,       , i8vals, INTEGER*8,        i8vals(*))
VARNM(put, IN,  text,   _all, text,   CHARACTER(len=*), text)
VARNM(put, IN,  int1,   _all, i1vals, INTEGER*1,        i1vals(*))
VARNM(put, IN,  int2,   _all, i2vals, INTEGER*2,        i2vals(*))
VARNM(put, IN,  int,    _all, ivals,  INTEGER,          ivals(*))
VARNM(put, IN,  real,   _all, rvals,  REAL,             rvals(*))
VARNM(put, IN,  double, _all, dvals,  DOUBLE PRECISION, dvals(*))
VARNM(put, IN,  int8,   _all, i8vals, INTEGER*8,        i8vals(*))
!
! End of varn subroutines
!

