dnl
dnl This utilities file contains common m4 macros used by C/Fortran program
dnl
dnl
define(`_CAT', `$1$2')dnl  concatenate two strings
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
define(`MPI2ITYPE',  `ifelse(`$1', `MPI_CHAR',               `text',
                             `$1', `MPI_SIGNED_CHAR',        `schar',
                             `$1', `MPI_UNSIGNED_CHAR',      `uchar',
                             `$1', `MPI_SHORT',              `short',
                             `$1', `MPI_UNSIGNED_SHORT',     `ushort',
                             `$1', `MPI_INT',                `int',
                             `$1', `MPI_UNSIGNED',           `uint',
                             `$1', `MPI_LONG',               `long',
                             `$1', `MPI_FLOAT',              `float',
                             `$1', `MPI_DOUBLE',             `double',
                             `$1', `MPI_LONG_LONG_INT',      `longlong',
                             `$1', `MPI_UNSIGNED_LONG_LONG', `ulonglong')')dnl
dnl
dnl dnl dnl
dnl
define(`ITYPE_LIST', `text, schar, uchar, short, ushort, int, uint, long, float, double, longlong, ulonglong')dnl
dnl
dnl
dnl dnl dnl
dnl
define(`CDF2_ITYPE_LIST', `text, schar, short, int, long, float, double')dnl
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
define(`NC_TYPE',`ifelse(
`$1', `text',      `NC_CHAR',
`$1', `schar',     `NC_BYTE',
`$1', `uchar',     `NC_UBYTE',
`$1', `short',     `NC_SHORT',
`$1', `ushort',    `NC_USHORT',
`$1', `int',       `NC_INT',
`$1', `long',      `NC_LONG',
`$1', `uint',      `NC_UINT',
`$1', `float',     `NC_FLOAT',
`$1', `double',    `NC_DOUBLE',
`$1', `longlong',  `NC_INT64',
`$1', `ulonglong', `NC_UINT64')')dnl
dnl
define(`NC_FILL_VALUE',`ifelse(
`$1', `text',      `NC_FILL_CHAR',
`$1', `schar',     `NC_FILL_BYTE',
`$1', `uchar',     `NC_FILL_UBYTE',
`$1', `short',     `NC_FILL_SHORT',
`$1', `ushort',    `NC_FILL_USHORT',
`$1', `int',       `NC_FILL_INT',
`$1', `long',      `NC_FILL_INT',
`$1', `uint',      `NC_FILL_UINT',
`$1', `float',     `NC_FILL_FLOAT',
`$1', `double',    `NC_FILL_DOUBLE',
`$1', `longlong',  `NC_FILL_INT64',
`$1', `ulonglong', `NC_FILL_UINT64')')dnl
dnl
define(`IFMT',`ifelse(
`$1', `text',      `%c',
`$1', `schar',     `%hhd',
`$1', `uchar',     `%hhu',
`$1', `short',     `%hd',
`$1', `ushort',    `%hu',
`$1', `int',       `%d',
`$1', `long',      `%ld',
`$1', `uint',      `%u',
`$1', `float',     `%g',
`$1', `double',    `%g',
`$1', `longlong',  `%lld',
`$1', `ulonglong', `%llu')')dnl
dnl
define(`PUT_VAR',`ifelse(
`$1', `text',      `ncmpi_put_var_$1_all',dnl
`$1', `schar',     `ncmpi_put_var_$1_all',dnl
`$1', `uchar',     `ncmpi_put_var_$1_all',dnl
`$1', `short',     `ncmpi_put_var_$1_all',dnl
`$1', `ushort',    `ncmpi_put_var_$1_all',dnl
`$1', `int',       `ncmpi_put_var_$1_all',dnl
`$1', `long',      `ncmpi_put_var_$1_all',dnl
`$1', `uint',      `ncmpi_put_var_$1_all',dnl
`$1', `float',     `ncmpi_put_var_$1_all',dnl
`$1', `double',    `ncmpi_put_var_$1_all',dnl
`$1', `longlong',  `ncmpi_put_var_$1_all',dnl
`$1', `ulonglong', `ncmpi_put_var_$1_all')')dnl
dnl
define(`GET_VAR',`ifelse(
`$1', `text',      `ncmpi_get_var_$1_all',dnl
`$1', `schar',     `ncmpi_get_var_$1_all',dnl
`$1', `uchar',     `ncmpi_get_var_$1_all',dnl
`$1', `short',     `ncmpi_get_var_$1_all',dnl
`$1', `ushort',    `ncmpi_get_var_$1_all',dnl
`$1', `int',       `ncmpi_get_var_$1_all',dnl
`$1', `long',      `ncmpi_get_var_$1_all',dnl
`$1', `uint',      `ncmpi_get_var_$1_all',dnl
`$1', `float',     `ncmpi_get_var_$1_all',dnl
`$1', `double',    `ncmpi_get_var_$1_all',dnl
`$1', `longlong',  `ncmpi_get_var_$1_all',dnl
`$1', `ulonglong', `ncmpi_get_var_$1_all')')dnl
dnl
define(`XTYPE_MAX',`ifelse(
`$1', `text',      `127',
`$1', `schar',     `127',
`$1', `uchar',     `255',
`$1', `short',     `32767',
`$1', `ushort',    `65535U',
`$1', `int',       `2147483647',
`$1', `long',      `2147483647',
`$1', `uint',      `4294967295U',
`$1', `float',     `3.402823466e+38f',
`$1', `double',    `1.79769313486230e+308',
`$1', `longlong',  `9223372036854775807LL',
`$1', `ulonglong', `18446744073709551615ULL')')dnl
dnl
