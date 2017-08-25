dnl
dnl This utilities file contains common m4 macros used by C/Fortran program
dnl
dnl
define(`_CAT', `$1$2')dnl  concatenate two strings
dnl
dnl
define(`NULL_CHAR', changequote([,])[changequote([,])'\0'changequote(`,')]changequote(`,'))
dnl
dnl dnl dnl
dnl
define(`NC2ITYPE', `ifelse(`$1', `text',      `char',
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
define(`SIZEOFITYPE', `ifelse(`$1', `text',      `SIZEOF_CHAR',
                              `$1', `schar',     `SIZEOF_SIGNED_CHAR',
                              `$1', `uchar',     `SIZEOF_UNSIGNED_CHAR',
                              `$1', `short',     `SIZEOF_SHORT',
                              `$1', `ushort',    `SIZEOF_UNSIGNED_SHORT',
                              `$1', `int',       `SIZEOF_INT',
                              `$1', `uint',      `SIZEOF_UNSIGNED_INT',
                              `$1', `long',      `SIZEOF_LONG',
                              `$1', `float',     `SIZEOF_FLOAT',
                              `$1', `double',    `SIZEOF_DOUBLE',
                              `$1', `longlong',  `SIZEOF_LONG_LONG',
                              `$1', `ulonglong', `SIZEOF_UNSIGNED_LONG_LONG')')dnl
dnl
dnl size of external NC data type is defined in CDF format specifications
dnl
define(`SIZEOFXTYPE', `ifelse(`$1', `NC_CHAR',   `1',
                              `$1', `NC_BYTE',   `1',
                              `$1', `NC_UBYTE',  `1',
                              `$1', `NC_SHORT',  `2',
                              `$1', `NC_USHORT', `2',
                              `$1', `NC_INT',    `4',
                              `$1', `NC_UINT',   `4',
                              `$1', `NC_LONG',   `4',
                              `$1', `NC_FLOAT',  `4',
                              `$1', `NC_DOUBLE', `8',
                              `$1', `NC_INT64',  `8',
                              `$1', `NC_UINT64', `8')')dnl
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
                             `$1', `ulonglong', `MPI_UNSIGNED_LONG_LONG',
                             `MPI_DATATYPE_NULL')')dnl
dnl
dnl dnl dnl
dnl
define(`NCTYPE2MPI', `ifelse(`$1', `NC_CHAR',   `MPI_CHAR',
                             `$1', `NC_BYTE',   `MPI_SIGNED_CHAR',
                             `$1', `NC_UBYTE',  `MPI_UNSIGNED_CHAR',
                             `$1', `NC_SHORT',  `MPI_SHORT',
                             `$1', `NC_USHORT', `MPI_UNSIGNED_SHORT',
                             `$1', `NC_INT',    `MPI_INT',
                             `$1', `NC_UINT',   `MPI_UNSIGNED',
                             `$1', `NC_LONG',   `MPI_LONG',
                             `$1', `NC_FLOAT',  `MPI_FLOAT',
                             `$1', `NC_DOUBLE', `MPI_DOUBLE',
                             `$1', `NC_INT64',  `MPI_LONG_LONG_INT',
                             `$1', `NC_UINT64', `MPI_UNSIGNED_LONG_LONG',
                             `MPI_DATATYPE_NULL')')dnl
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
define(`CollIndep', `ifelse(`$1', `_all', `NC_REQ_COLL', `NC_REQ_INDEP')')dnl
define(`ReadWrite', `ifelse(`$1',  `get', `NC_REQ_WR',
                            `$1', `iget', `NC_REQ_RD',
                                          `NC_REQ_WR')')dnl
define(`BufArgs',   `ifelse(`$2', `',
                            `ifelse($1, `get', `void *buf,', `const void *buf,')
                             MPI_Offset bufcount, MPI_Datatype buftype',
                            `ifelse($1, `get',       `NC2ITYPE($2) *buf',
                                               `const NC2ITYPE($2) *buf')')')
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
define(`ArgStartCountStrideMap', `ifelse(
       `$1', `',  `NULL,  NULL,  NULL,   NULL',
       `$1', `1', `start, NULL,  NULL,   NULL',
       `$1', `a', `start, count, NULL,   NULL',
       `$1', `s', `start, count, stride, NULL',
       `$1', `m', `start, count, stride, imap')')dnl
dnl
define(`ArgStartCountStride', `ifelse(
       `$1', `',  `NULL,  NULL,  NULL',
       `$1', `1', `start, NULL,  NULL',
       `$1', `a', `start, count, NULL',
                  `start, count, stride')')dnl
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
define(`PUT_VAR',`ifdef(`PNETCDF',`ncmpi_put_var_$1_all',`nc_put_var_$1')')dnl
dnl
define(`GET_VAR',`ifdef(`PNETCDF',`ncmpi_get_var_$1_all',`nc_get_var_$1')')dnl
dnl
define(`PUT_VAR1',`ifdef(`PNETCDF',`ncmpi_put_var1_$1_all',`nc_put_var1_$1')')dnl
dnl
define(`GET_VAR1',`ifdef(`PNETCDF',`ncmpi_get_var1_$1_all',`nc_get_var1_$1')')dnl
dnl
define(`PUT_VARA',`ifdef(`PNETCDF',`ncmpi_put_vara_$1_all',`nc_put_vara_$1')')dnl
dnl
define(`GET_VARA',`ifdef(`PNETCDF',`ncmpi_get_vara_$1_all',`nc_get_vara_$1')')dnl
dnl
define(`PUT_VARS',`ifdef(`PNETCDF',`ncmpi_put_vars_$1_all',`nc_put_vars_$1')')dnl
dnl
define(`GET_VARS',`ifdef(`PNETCDF',`ncmpi_get_vars_$1_all',`nc_get_vars_$1')')dnl
dnl
define(`PUT_VARM',`ifdef(`PNETCDF',`ncmpi_put_varm_$1_all',`nc_put_varm_$1')')dnl
dnl
define(`GET_VARM',`ifdef(`PNETCDF',`ncmpi_get_varm_$1_all',`nc_get_varm_$1')')dnl
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
