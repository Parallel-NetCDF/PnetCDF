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
