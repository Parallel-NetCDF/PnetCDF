!
! Begin of varn subroutines:
!
! Note that we use starts(*) and counts(*) instead of starts(:,:) and
! counts(:,:), because we will rearrange the Fortran order for the
! dimensions to C order and would like the arrays in a contiguous space
!
! $Id$
!

define(`VARN_FLEXIBLE',dnl
`dnl
    INTEGER FUNCTION nfmpi_$1_varn$3(ncid, varid, num, starts, counts, buf, bufcount, buftype)
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     <type>,                        INTENT($2)  :: buf(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
                     INTEGER,                       INTENT(IN)  :: buftype
    END FUNCTION nfmpi_$1_varn$3
')dnl
#if 0
!
! flexible APIs, not ready yet for Fortran90
!
VARN_FLEXIBLE(get, OUT)
VARN_FLEXIBLE(get, OUT, _all)
VARN_FLEXIBLE(put, IN)
VARN_FLEXIBLE(put, IN,  _all)
#endif

define(`VARN',dnl
`dnl
    INTEGER FUNCTION nfmpi_$1_varn_$3$4(ncid, varid, num, starts, counts, $5)
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     $6,			    INTENT($2)  :: $7
    END FUNCTION nfmpi_$1_varn_$3$4
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
!
! End of varn subroutines:
!

