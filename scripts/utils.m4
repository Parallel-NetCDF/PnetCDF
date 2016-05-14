dnl
dnl This utilities file contains common m4 macros used by C/Fortran program
dnl
dnl
dnl dnl dnl
dnl
define(`FUNC2ITYPE', `ifelse(`$1', `text',      `char',
                             `$1', `schar',     `schar',
                             `$1', `uchar',     `uchar',
                             `$1', `short',     `short',
                             `$1', `ushort',    `ushort',
                             `$1', `int',       `int',
                             `$1', `uint',      `uint',
                             `$1', `long',      `long',
                             `$1', `float',     `float',
                             `$1', `double',    `double',
                             `$1', `longlong',  `long long',
                             `$1', `ulonglong', `unsigned long long')')dnl
dnl
dnl dnl dnl
dnl
define(`ITYPE2MPI',  `ifelse(`$1', `text',      `MPI_CHAR',
                             `$1', `schar',     `MPI_SIGNED_CHAR',
                             `$1', `uchar',     `MPI_UNSIGNED_CHAR',
                             `$1', `short',     `MPI_SHORT',
                             `$1', `ushort',    `MPI_UNSIGNED_SHORT',
                             `$1', `int',       `MPI_INT',
                             `$1', `uint',      `MPI_UNSIGNED',
                             `$1', `long',      `MPI_LONG',
                             `$1', `float',     `MPI_FLOAT',
                             `$1', `double',    `MPI_DOUBLE',
                             `$1', `longlong',  `MPI_LONG_LONG_INT',
                             `$1', `ulonglong', `MPI_UNSIGNED_LONG_LONG')')dnl
dnl
dnl dnl dnl
dnl
define(`ITYPE_LIST', `text, schar, uchar, short, ushort, int, uint, long, float, double, longlong, ulonglong')dnl
dnl
dnl
define(`CollIndep', `ifelse(`$1', `_all', `COLL_IO', `INDEP_IO')')dnl
define(`ReadWrite', `ifelse(`$1',  `get', `READ_REQ',
                            `$1', `iget', `READ_REQ',
                                          `WRITE_REQ')')dnl
define(`BufArgs',   `ifelse(`$2', `',
                            `ifelse($1, `put', `const void *buf,', `void *buf,')
                             MPI_Offset   bufcount,
                             MPI_Datatype buftype',
                            `ifelse($1, `put', `const FUNC2ITYPE($2) *buf',
                                                     `FUNC2ITYPE($2) *buf')')')
dnl
dnl index arguments for APIs of different kinds
dnl
define(`ArgKind', `ifelse(
       `$1', `1', `const MPI_Offset *start,',
       `$1', `a', `const MPI_Offset *start,
                   const MPI_Offset *count,',
       `$1', `s', `const MPI_Offset *start,
                   const MPI_Offset *count,
                   const MPI_Offset *stride,',
       `$1', `m', `const MPI_Offset *start,
                   const MPI_Offset *count,
                   const MPI_Offset *stride,
                   const MPI_Offset *imap,')')dnl
dnl
dnl arguments passed to a function for APIs of different kinds
dnl
define(`ArgStartCount', `ifelse(
       `$1', `',  `NULL,  NULL',
       `$1', `1', `start, NULL',
                  `start, count')')dnl
dnl
define(`ArgStrideMap', `ifelse(
       `$1', `s', `stride, NULL',
       `$1', `m', `stride, imap',
                  `NULL, NULL')')dnl
dnl
define(`API_KIND', `ifelse(
       `$1', `1', `API_VAR1',
       `$1', `a', `API_VARA',
       `$1', `s', `API_VARS',
       `$1', `m', `API_VARM',
       `$1', `d', `API_VARD',
       `$1', `n', `API_VARN',
       `$1', `',  `API_VAR')')dnl
dnl
