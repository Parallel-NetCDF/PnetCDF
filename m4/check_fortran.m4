dnl Check for Fortran-90 compiler.
dnl
AC_DEFUN([UD_PROG_F90],
[
    case "${F90+set}" in
	set)
	    AC_MSG_CHECKING(user-defined Fortran-90 compiler \"$F90\")
            AC_LANG_PUSH([Fortran])
            AC_COMPILE_IFELSE(
               [AC_LANG_SOURCE([
		   subroutine foo(bar)
		   integer, intent(in) :: bar
		   end subroutine foo
               ])],
               [AC_MSG_RESULT(works)],
               [AC_MSG_RESULT(failed to compile test program)
		unset F90
               ]
            )
            AC_LANG_POP([Fortran])
	    ;;
	*)
	    case "${FC+set}" in
		set)
		    F90=$FC
		    F90FLAGS="${F90FLAGS-${FFLAGS--O}}"
		    F90LIBS="${F90LIBS-${FLIBS}}"
		    cat <<EOF >conftest.f90
			program foo
			call bar(1)
			end
			subroutine bar(bof)
			integer, intent(in) :: bof
			end
EOF
		    AC_MSG_CHECKING(\"$F90\" as Fortran-90 compiler)
		    doit='$F90 -o conftest ${F90FLAGS} conftest.f90 ${F90LIBS}'
		    if AC_TRY_EVAL(doit); then
			doit=./conftest
			if AC_TRY_EVAL(doit); then
			    AC_MSG_RESULT(works)
			else
			    AC_MSG_RESULT(failed to build executable program)
			    unset F90
			fi
		    else
			AC_MSG_RESULT(failed to build test program)
			unset F90
		    fi
		    ${RM} -f conftest*
		    ;;
	    esac
	    case "${F90-unset}" in
		unset)
		    cat <<EOF >conftest.f90
			program foo
			call bar(1)
			end
			subroutine bar(bof)
			integer, intent(in) :: bof
			end
EOF
		    for f90 in xlf90 pgf90 f90; do
			AC_CHECK_PROG(F90, $f90, $f90)
			case "${F90}" in
			    '')
				;;
			    *)
				AC_MSG_CHECKING(Fortran-90 compiler \"$F90\")
				doit='$F90 -o conftest ${F90FLAGS} conftest.f90 ${F90LIBS}'
				if AC_TRY_EVAL(doit); then
				    doit=./conftest
				    if AC_TRY_EVAL(doit); then
					AC_MSG_RESULT(works)
					break
				    else
					AC_MSG_RESULT(
					    failed to build executable program)
					unset F90
					unset ac_cv_prog_F90
				    fi
				else
				    AC_MSG_RESULT(failed to build test program)
				    unset F90
				    unset ac_cv_prog_F90
				fi
				;;
			esac
		    done
		    ${RM} -f conftest*
		    case "${F90}" in
			'') AC_MSG_WARN(
			    "Could not find working Fortran-90 compiler")
			    ;;
		    esac
		    ;;
	    esac
	    ;;
    esac
    case "${F90}" in
	'')
	    AC_MSG_WARN("The Fortran-90 interface will not be built")
	    ;;
	*f90*)
	    case `uname -s` in
		IRIX*)
		    NETCDF_MOD=NETCDF.mod
		    ;;
		*)
		    NETCDF_MOD=netcdf.mod
		    ;;
	    esac
	    AC_SUBST(NETCDF_MOD)
	    ;;
    esac
    AC_SUBST(F90)
    AC_SUBST(F90FLAGS)
    AC_SUBST(F90LIBS)
])

dnl Check for Fortran-77 compiler.
dnl
AC_DEFUN([UD_PROG_FC],
[
    AC_BEFORE([UD_FORTRAN_TYPES])
    case "${FC+set}" in
	set)
	    case "$FC" in
		'')
		    AC_MSG_WARN(Fortran-77 compiler is explicitly set to null)
		    ;;
		*)
		    AC_MSG_CHECKING(user-defined Fortran-77 compiler \"$FC\")
                    AC_LANG_PUSH([Fortran 77])
                    AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                       [AC_MSG_RESULT(works)],
                       [AC_MSG_RESULT(failed to compile test program)
                        FC=
                       ]
                    )
                    AC_LANG_POP([Fortran 77])
		    ;;
	    esac
	    ;;
	*)
            dnl FC is not set, let's try F90 as FC if F90 is set
	    case "${F90+set}" in
		set)
		    FC=$F90
		    FFLAGS="${FFLAGS-${F90FLAGS--O}}"
		    FLIBS="${FLIBS-${F90LIBS-}}"
		    AC_MSG_CHECKING(\"$FC\" as Fortran-77 compiler)
                    AC_LANG_PUSH([Fortran])
                    AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                       [AC_MSG_RESULT(works)],
                       [AC_MSG_RESULT(failed to compile test program)
			unset FC
                       ]
                    )
                    AC_LANG_POP([Fortran])
		    ;;
		*)
		    ;;
	    esac
	    case "${FC-unset}" in
		unset)
		    case `uname -sr` in
			AIX*)
			    # xlf90(1) thinks fortran/ftest.F has bad syntax.
			    forts="xlf f77 gfortran"
			    ;;
			BSD/OS*|FreeBSD*)
			    forts="f77 fort77 g77 gfortran"
			    ;;
			HP-UX*)
			    # f77(1) doesn't have the -L option.
			    forts=fort77
			    FLIBS=-lU77
			    ;;
			IRIX*)
			    # f90(1) can't link with c89(1)-compiled objects
			    forts=f77
			    ;;
			IRIX64*)
			    forts='f77 g77 gfortran fort77 '
			    ;;
			Linux*)
			    forts="pgf90 f77 fort77 g77 gfortran"
			    ;;
			OSF1*)
			    # The use of f90(1) results in the following for
			    # an unknown reason (`make' works in the fortran/
			    # directory):
			    # f90 -c -I../libsrc ftest.F
			    # Last chance handler: pc = 0xa971b8,
			    # sp = 0x3fece0, ra = 0xa971b8
			    # Last chance handler: internal exception: unwinding
			    forts="f77 gfortran"
			    ;;
			'SunOS 4'*)
			    forts='f77 g77 gfortran fort77'
			    ;;
			'SunOS 5'*)
			    # SunOS's f90(1) has problems passing a C `char'
			    # as a Fortran `integer*1' => use f77(1)
			    forts="pgf90 f77"
			    ;;
			sn*|UNICOS*|unicos*)
			    forts="fort77 cf77 f77 g77 gfortran f90"
			    ;;
			*)
			    forts="xlf fort77 ghf77 f77 cf77 g77 gfortran xlf90 f90"
			    ;;
		    esac
		    for fc in $forts; do
			AC_CHECK_PROG(FC, $fc, $fc)
			case "${FC}" in
			    '')
				;;
			    *)
				#
				# On some systems, a discovered compiler
				# nevertheless won't work (due to licensing,
				# for example); thus, we check the compiler
				# with a test program.
				#
                                AC_LANG_PUSH([Fortran])
                                AC_COMPILE_IFELSE([AC_LANG_CALL([],[FOO])],
                                   [AC_MSG_RESULT(works)],
                                   [AC_MSG_RESULT(failed to compile test program)
			            unset FC
				    unset ac_cv_prog_FC
                                   ]
                                )
                                AC_LANG_POP([Fortran])
				;;
			esac
		    done
		    ${RM} -f conftest.*
		    case "${FC}" in
			'') AC_MSG_WARN(
				"Could not find working Fortran-77 compiler")
			    ;;
		    esac
		    ;;
	    esac
	    ;;
    esac
    case "${FC}" in
	'') AC_MSG_WARN("The Fortran-77 interface will not be built")
	    ;;
    esac
    AC_SUBST(FC)
    AC_SUBST(FFLAGS)
    AC_SUBST(FLIBS)
    #
    # Set the make(1) macro for compiling a .F file.
    #
    case "${FPP-}" in
    '')
	AC_MSG_CHECKING(for Fortran .F compiler)
	AC_MSG_RESULT($COMPILE_F)
	case "${COMPILE_F-unset}" in
	unset)
	    case "${FC}" in
	    '')
		COMPILE_F=
		;;
	    *)
		cat >conftest.h <<\EOF
#define J 1
EOF
                AC_LANG_PUSH([Fortran])
                AC_FC_SRCEXT([F])
		AC_MSG_CHECKING(if Fortran-77 compiler handles *.F files)
                AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#include "conftest.h"
#define N 5
		  real r(J,N)
                   ])],
                   [AC_MSG_RESULT(yes)
                    COMPILE_F='$(COMPILE.f)'],
                   [AC_MSG_RESULT(no)
		    COMPILE_F=
                   ]
                )
                AC_FC_SRCEXT([f])
                AC_LANG_POP([Fortran])
		${RM} -f conftest.h
		;;
	    esac
	    ;;
	esac
	;;
    *)
	unset COMPILE_F
	;;
    esac
    case "${COMPILE_F-}" in
	'') UD_PROG_FPP;;
    esac
    AC_SUBST(COMPILE_F)
    FPPFLAGS=${FPPFLAGS-}
    AC_SUBST(FPPFLAGS)
])

dnl Check for Fortran preprocessor.
dnl
AC_DEFUN([UD_PROG_FPP],
[
    AC_MSG_CHECKING(for Fortran preprocessor)
    case "$FPP" in
    '')
	AC_REQUIRE([AC_PROG_CPP])
	FPP="$CPP"
	;;
    esac
    AC_MSG_RESULT($FPP)
    AC_SUBST(FPP)
])

dnl Check for a Fortran type equivalent to a netCDF type.
dnl
dnl UD_CHECK_FORTRAN_NCTYPE(forttype, possibs, nctype)
dnl
AC_DEFUN([UD_CHECK_FORTRAN_NCTYPE],
[
    AC_MSG_CHECKING([for Fortran-equivalent to netCDF "$3"])
dnl     for type in $2; do
dnl         cat >conftest.f <<EOF
dnl                $type foo
dnl                end
dnl EOF
dnl         doit='$FC -c ${FFLAGS} conftest.f'
dnl         if AC_TRY_EVAL(doit); then
dnl             break
dnl         fi
dnl     done
dnl     ${RM} -f conftest.f conftest.o

    AC_LANG_PUSH([Fortran])
    for type in $2; do
        AC_COMPILE_IFELSE(
           [AC_LANG_SOURCE([
               $type foo
               end
           ])],
           [break]
        )
    done
    AC_LANG_POP([Fortran])
    AC_DEFINE_UNQUOTED($1, $type)
    AC_MSG_RESULT($type)
    $1=$type
])

dnl Check for a Fortran type equivalent to a C type.
dnl
dnl UD_CHECK_FORTRAN_CTYPE(v3forttype, v2forttype, ctype, min, max)
dnl
AC_DEFUN([UD_CHECK_FORTRAN_CTYPE],
[
    AC_MSG_CHECKING([for Fortran-equivalent to C "$3"])
    AC_LANG_PUSH([Fortran])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
        subroutine sub(values, minval, maxval)
        implicit        none
        $2              values(5), minval, maxval
        minval = values(2)
        maxval = values(4)
        if (values(2) .ge. values(4)) then
            minval = values(4)
            maxval = values(2)
        endif
        end
       ])],
       [ac_found_f2c=yes], [ac_found_f2c=no]
    )
    AC_LANG_POP([Fortran])
    if test "x${ac_found_f2c}" = xyes ; then
dnl     cat >conftest.f <<EOF
dnl         subroutine sub(values, minval, maxval)
dnl         implicit        none
dnl         $2              values(5), minval, maxval
dnl         minval = values(2)
dnl         maxval = values(4)
dnl         if (values(2) .ge. values(4)) then
dnl             minval = values(4)
dnl             maxval = values(2)
dnl         endif
dnl         end
dnl EOF
dnl     doit='$FC -c ${FFLAGS} conftest.f'
dnl     if AC_TRY_EVAL(doit); then
dnl         mv conftest.o conftestf.o
	cat >conftest.c <<EOF
#include <limits.h>
#include <float.h>
void main()
{
$3		values[[]] = {0, $4, 0, $5, 0};
$3		minval, maxval;
int	$FCALLSCSUB($3*, $3*, $3*);
$FCALLSCSUB(values, &minval, &maxval);
return(!(minval == $4 && maxval == $5));
}
EOF
	doit='$CC -o conftest ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} conftest.c conftestf.o ${LIBS}'
	if AC_TRY_EVAL(doit); then
	    doit=./conftest
	    if AC_TRY_EVAL(doit); then
		AC_MSG_RESULT($2)
		$1=$2
		AC_DEFINE_UNQUOTED($1,$2)
	    else
		AC_MSG_RESULT(no equivalent type)
		unset $1
	    fi
	else
	    AC_MSG_ERROR(Could not compile-and-link conftest.c and conftestf.o)
	fi
    else
	AC_MSG_ERROR(Could not compile conftest.f)
    fi
    ${RM} -f conftest*
    unset ac_found_f2c
])

dnl Check for a Fortran data type.
dnl
dnl UD_CHECK_FORTRAN_TYPE(varname, ftypes)
dnl
AC_DEFUN([UD_CHECK_FORTRAN_TYPE],
[
    AC_LANG_PUSH([Fortran])
    for ftype in $2; do
	AC_MSG_CHECKING([for Fortran "$ftype"])
        AC_COMPILE_IFELSE(
           [AC_LANG_SOURCE([
               subroutine sub(value)
               $ftype value
               end
           ])],
           [AC_MSG_RESULT(yes)
	    $1=$ftype
	    AC_DEFINE_UNQUOTED([$1], [$ftype])
            break],
           [AC_MSG_RESULT(no)]
        )
    done
    AC_LANG_POP([Fortran])
])

dnl Check for the name format of a Fortran-callable C routine.
dnl
dnl UD_CHECK_FCALLSCSUB
AC_DEFUN([UD_CHECK_FCALLSCSUB],
[
    dnl AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
	'') ;;
	*)
	    AC_REQUIRE([UD_PROG_NM])
	    AC_BEFORE([UD_CHECK_FORTRAN_CTYPE])
	    AC_BEFORE([UD_CHECK_CTYPE_FORTRAN])
	    AC_MSG_CHECKING([for C-equivalent to Fortran routine "SUB"])
            AC_FC_FUNC([sub], [FCALLSCSUB])
            AC_MSG_RESULT([$FCALLSCSUB])
	    ;;
    esac
])

dnl Check for a C type equivalent to a Fortran type.
dnl
dnl UD_CHECK_CTYPE_FORTRAN(ftype, ctypes, fmacro_root)
dnl
AC_DEFUN([UD_CHECK_CTYPE_FORTRAN],
[
    cat >conftestf.f <<EOF
           $1 values(4)
           integer status, sub
           data values /-1, -2, -3, -4/
           status = sub(values)
           end
EOF
    ac_cv_ctype_fortran=no
    AC_MSG_CHECKING([if Fortran "$1" is ])
    for ctype in $2; do
	dnl AC_MSG_CHECKING(if Fortran \"$1\" is C \"$ctype\")
	cat >conftest.c <<EOF
	    int $FCALLSCSUB($ctype values[[4]])
	    {
		return(values[[1]] != -2 || values[[2]] != -3);
	    }
EOF
	doit='$CC -c ${CPPFLAGS} ${CFLAGS} conftest.c'
	if AC_TRY_EVAL(doit); then
	    doit='$FC ${FFLAGS} -c conftestf.f'
	    if AC_TRY_EVAL(doit); then
	        doit='$FC -o conftest ${FFLAGS} ${FLDFLAGS} conftestf.o conftest.o ${LDFLAGS} ${LIBS}'
	        if AC_TRY_EVAL(doit); then
		    doit=./conftest
		    if AC_TRY_EVAL(doit); then
		        dnl AC_MSG_RESULT(yes)
		        AC_MSG_RESULT(["$ctype" in C])
		        cname=`echo $ctype | tr ' abcdefghijklmnopqrstuvwxyz' \
			    _ABCDEFGHIJKLMNOPQRSTUVWXYZ`
		        AC_DEFINE_UNQUOTED(NF_$3[]_IS_C_$cname)
                        ac_cv_ctype_fortran=yes
		        break
		    fi
	        else
		    AC_MSG_ERROR([Could not link conftestf.o and conftest.o])
	        fi
	    else
		AC_MSG_ERROR([Could not compile conftestf.f])
	    fi
	else
	    AC_MSG_ERROR([Could not compile conftest.c])
	fi
    done
    ${RM} -f conftest*

    if test "$ac_cv_ctype_fortran" = no ; then
        AC_MSG_RESULT(no correspond data type in C)
    fi
    unset ac_cv_ctype_fortran
])

dnl Get information about Fortran data types.
dnl
AC_DEFUN([UD_FORTRAN_TYPES],
[
    dnl AC_REQUIRE([UD_PROG_FC])
    case "$FC" in
    '')
	;;
    *)
	AC_REQUIRE([UD_CHECK_FCALLSCSUB])
	UD_CHECK_FORTRAN_TYPE([NF_INT1_T], [integer*1 byte "integer(kind=1)"])
	UD_CHECK_FORTRAN_TYPE([NF_INT2_T], [integer*2 "integer(kind=2)"])
	UD_CHECK_FORTRAN_TYPE([NF_INT8_T], [integer*8 "integer(kind=8)"])

	case "${NF_INT1_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT1_T, "signed char" short int long, INT1)
		;;
	esac
	case "${NF_INT2_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT2_T, short int long, INT2)
		;;
	esac
	case "${NF_INT8_T}" in
	    '') ;;
	    *)  UD_CHECK_CTYPE_FORTRAN($NF_INT8_T, int long "long long", INT8)
		;;
	esac
	UD_CHECK_CTYPE_FORTRAN(integer, int long, INT)
	UD_CHECK_CTYPE_FORTRAN(real, float double, REAL)
	UD_CHECK_CTYPE_FORTRAN(doubleprecision, double float, DOUBLEPRECISION)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCBYTE_T, byte integer*1 integer, byte)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCSHORT_T, integer*2 integer, short)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_SHORT_T, $NCSHORT_T, short, SHRT_MIN, SHRT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCLONG_T, integer*4 integer, long)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_INT_T, integer, int, INT_MIN, INT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCFLOAT_T, real*4 real, float)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_FLOAT_T, $NCFLOAT_T, float, FLT_MIN, FLT_MAX)

dnl	UD_CHECK_FORTRAN_NCTYPE(NCDOUBLE_T, real*8 doubleprecision real, double)
dnl	UD_CHECK_FORTRAN_CTYPE(NF_DOUBLE_T, $NCDOUBLE_T, double, DBL_MIN, DBL_MAX)
	;;
    esac
])

dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_SRCEXT],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile preprocessed .$1 files],
                ac_cv_fc_pp_srcext_$1,
[ac_ext=$1
FCFLAGS_SRCEXT_save=$FCFLAGS_SRCEXT
FCFLAGS_SRCEXT=
ac_fcflags_pp_srcext_save=$ac_fcflags_srcext
ac_fcflags_srcext=
ac_cv_fc_pp_srcext_$1=unknown
case $ac_ext in #(
  [[fF]]77) ac_try=f77-cpp-input;; #(
  [[fF]]) ac_try=f77-cpp-input;; #(
  *) ac_try=f95-cpp-input;;
esac
for ac_flag in none -ftpp -fpp -Tf "-fpp -Tf" -xpp=fpp -Mpreprocess "-e Z" \
               -cpp -xpp=cpp -qsuffix=cpp=$1 "-x $ac_try" +cpp -Cpp; do
  test "x$ac_flag" != xnone && FCFLAGS_SRCEXT="$ac_flag" && ac_fcflags_srcext="$ac_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 0
#include <ac_nonexistent.h>
      choke me
#endif]])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 1
#include <ac_nonexistent.h>
      choke me
#endif]])],
       [],
       [ac_cv_fc_pp_srcext_$1=$ac_flag; break])])
done
${RM} -f conftest.$ac_objext conftest.$1
FCFLAGS_SRCEXT=$FCFLAGS_SRCEXT_save
ac_fcflags_srcext=$ac_fcflags_pp_srcext_save
])
if test "x$ac_cv_fc_pp_srcext_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran could not compile preprocessed .$1 files])])
else
  ac_fc_srcext=$1
  if test "x$ac_cv_fc_pp_srcext_$1" = xnone; then
    ac_fcflags_srcext=""
    FCFLAGS_[]$1[]=""
  else
    ac_fcflags_srcext=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_[]$1[]=$ac_cv_fc_pp_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
AC_LANG_POP(Fortran)dnl
])# UD_FC_PP_SRCEXT

dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_DEFINE],
[AC_LANG_PUSH([Fortran])dnl
ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_ext_saved=$ac_ext
UD_FC_PP_SRCEXT([F])
ac_ext=F
AC_CACHE_CHECK([how to define symbols for preprocessed Fortran],
  [ac_cv_fc_pp_define],
[ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_cv_fc_pp_define=unknown
ac_fc_pp_define_FCFLAGS_save=$FCFLAGS
for ac_flag in -D -WF,-D -Wp,-D -Wc,-D
do
  FCFLAGS="$ac_fc_pp_define_FCFLAGS_save ${ac_flag}FOOBAR ${ac_flag}ZORK=42"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#ifndef FOOBAR
      choke me
#endif
#if ZORK != 42
      choke me
#endif]])],
    [ac_cv_fc_pp_define=$ac_flag])
  test x"$ac_cv_fc_pp_define" != xunknown && break
done
FCFLAGS=$ac_fc_pp_define_FCFLAGS_save
])
ac_fc_srcext=$ac_fc_pp_define_srcext_save
if test "x$ac_cv_fc_pp_define" = xunknown; then
  FC_DEFINE=
  m4_default([$2],
             [AC_MSG_ERROR([Fortran does not allow to define preprocessor symbols], 77)])
else
  FC_DEFINE=$ac_cv_fc_pp_define
  $1
fi
ac_ext=$ac_ext_saved
AC_SUBST([FC_DEFINE])dnl
AC_LANG_POP([Fortran])dnl
])

# AC_PROG_FC_MOD
# ---------------
dnl Note that Mac OSX file system is case-insensitive, so this function does
dnl not work precisely on Mac. Hence, we check whether the file system can
dnl find the file with name in lowercase. If not, we say the mod file is in
dnl uppercase.
AC_DEFUN([UD_PROG_FC_UPPERCASE_MOD],
[
AC_REQUIRE([UD_FC_MODULE_EXTENSION])
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether Fortran 90 compiler capitalizes .mod filenames])
AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([
        module conftest
        end module conftest
    ])]
)
dnl ac_try='$F90 ${F90FLAGS} conftest.f90 ${F90LIBS}>&AS_MESSAGE_LOG_FD'
dnl AC_TRY_EVAL(ac_try)
if test -f conftest.${FC_MODEXT} ; then
   ac_cv_prog_f90_uppercase_mod=no
else
   ac_cv_prog_f90_uppercase_mod=yes
   ${RM} -f CONFTEST.${FC_MODEXT}
fi
AC_MSG_RESULT($ac_cv_prog_f90_uppercase_mod)
${RM} -f conftest*
AC_LANG_POP(Fortran)
])

dnl steal from autoconf 2.69
# UD_FC_MODULE_EXTENSION
# ----------------------
# Find the Fortran 90 module file extension.  The module extension is stored
# in the variable FC_MODEXT and empty if it cannot be determined.  The result
# or "unknown" is cached in the cache variable ac_cv_fc_module_ext.
AC_DEFUN([UD_FC_MODULE_EXTENSION],
[AC_CACHE_CHECK([Fortran 90 module extension], [ac_cv_fc_module_ext],
[AC_LANG_PUSH(Fortran)
mkdir conftest.dir
cd conftest.dir
ac_cv_fc_module_ext=unknown
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [ac_cv_fc_module_ext=`ls | sed -n 's,conftest_module\.,,p'`
   if test x$ac_cv_fc_module_ext = x; then
dnl Some F90 compilers use upper case characters for the module file name.
     ac_cv_fc_module_ext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
   fi])
cd ..
${RM} -rf conftest.dir
AC_LANG_POP(Fortran)
])
FC_MODEXT=$ac_cv_fc_module_ext
if test "$FC_MODEXT" = unknown; then
  FC_MODEXT=
fi
AC_SUBST([FC_MODEXT])dnl
])

dnl steal from autoconf 2.69
dnl Fix a bug that mistakenly sets FC_MODINC to -M when Fujitsu frtpx is used.
dnl The correct one should be -I. The fix is shown below.
dnl bug<    for ac_flag in -M -I '-I ' '-M ' -p '-mod ' '-module ' '-Am -I'; do
dnl ---
dnl fix>    for ac_flag in -I '-I ' -M '-M ' -p '-mod ' '-mdir ' '-module ' '-Am -I'; do
dnl
# UD_FC_MODULE_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ---------------------------------------------------------------------
# Find a flag to include Fortran 90 modules from another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_flag variable.
# It may contain significant trailing whitespace.
#
# Known flags:
# gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated), -Jdir for writing)
# g95: -I dir (-fmod=dir for writing)
# SUN: -Mdir, -M dir (-moddir=dir for writing;
#                     -Idir for includes is also searched)
# HP: -Idir, -I dir (+moddir=dir for writing)
# IBM: -Idir (-qmoddir=dir for writing)
# Intel: -Idir -I dir (-mod dir for writing)
# Absoft: -pdir
# Lahey: -mod dir
# Cray: -module dir, -p dir (-J dir for writing)
#       -e m is needed to enable writing .mod files at all
# Compaq: -Idir
# NAGWare: -I dir (-mdir dir for writing)
# PathScale: -I dir  (but -module dir is looked at first)
# Portland: -module dir (first -module also names dir for writing)
# Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.', then -I)
#                    (-Am indicates how module information is saved)
AC_DEFUN([UD_FC_MODULE_FLAG],[
AC_CACHE_CHECK([Fortran 90 module inclusion flag], [ac_cv_fc_module_flag],
[AC_LANG_PUSH([Fortran])
ac_cv_fc_module_flag=unknown
mkdir conftest.dir
cd conftest.dir
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [cd ..
   ac_fc_module_flag_FCFLAGS_save=$FCFLAGS
   # Flag ordering is significant for gfortran and Sun.
   for ac_flag in -I '-I ' -M '-M ' -p '-mod ' '-mdir ' '-module ' '-Am -I'; do
     # Add the flag twice to prevent matching an output flag.
     FCFLAGS="$ac_fc_module_flag_FCFLAGS_save ${ac_flag}conftest.dir ${ac_flag}conftest.dir"
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ac_cv_fc_module_flag="$ac_flag"])
     if test "$ac_cv_fc_module_flag" != unknown; then
       break
     fi
   done
   FCFLAGS=$ac_fc_module_flag_FCFLAGS_save
])
${RM} -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_flag" != unknown; then
  FC_MODINC=$ac_cv_fc_module_flag
  $1
else
  FC_MODINC=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag for module search path])])
fi
AC_SUBST([FC_MODINC])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODINC in #(
  *\ ) FC_MODINC=$FC_MODINC'${ac_empty}' ;;
esac])dnl
])


dnl steal from autoconf 2.69
# UD_FC_MODULE_OUTPUT_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------------
# Find a flag to write Fortran 90 module information to another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_output_flag variable.
# It may contain significant trailing whitespace.
#
# For known flags, see the documentation of AC_FC_MODULE_FLAG above.
AC_DEFUN([UD_FC_MODULE_OUTPUT_FLAG],[
AC_CACHE_CHECK([Fortran 90 module output flag], [ac_cv_fc_module_output_flag],
[AC_LANG_PUSH([Fortran])
mkdir conftest.dir conftest.dir/sub
cd conftest.dir
ac_cv_fc_module_output_flag=unknown
ac_fc_module_output_flag_FCFLAGS_save=$FCFLAGS
# Flag ordering is significant: put flags late which some compilers use
# for the search path.
for ac_flag in -J '-J ' -fmod= -moddir= +moddir= -qmoddir= '-mod ' \
	      '-mdir ' '-module ' -M '-Am -M' '-e m -J '; do
  FCFLAGS="$ac_fc_module_output_flag_FCFLAGS_save ${ac_flag}sub"
  AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
    [cd sub
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end]],
       [ac_cv_fc_module_output_flag="$ac_flag"])
     cd ..
     if test "$ac_cv_fc_module_output_flag" != unknown; then
       break
     fi])
done
FCFLAGS=$ac_fc_module_output_flag_FCFLAGS_save
cd ..
${RM} -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_output_flag" != unknown; then
  FC_MODOUT=$ac_cv_fc_module_output_flag
  $1
else
  FC_MODOUT=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag to write module information to])])
fi
AC_SUBST([FC_MODOUT])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODOUT in #(
  *\ ) FC_MODOUT=$FC_MODOUT'${ac_empty}' ;;
esac])dnl
])

dnl
dnl steal from autoconf 2.69
dnl customized AC_FC_FREEFORM: to just get ac_cv_fc_freeform without
dnl appending ac_cv_fc_freeform to FCFLAGS
dnl
# AC_FC_FREEFORM([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept
# free-format source code, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using new extension) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
#
# The known flags are:
#        -ffree-form: GNU g77, gfortran, g95
#         -FR, -free: Intel compiler (icc, ecc, ifort)
#              -free: Compaq compiler (fort), Sun compiler (f95)
# -qfree=f90, -qfree: IBM compiler (xlf)
# -Mfree, -Mfreeform: Portland Group compiler
#          -freeform: SGI compiler
#        -8, -f free: Absoft Fortran
#       +source=free: HP Fortran
#    (-)-nfix, -Free: Lahey/Fujitsu Fortran
#              -free: NAGWare
#         -f, -Wf,-f: f2c (but only a weak form of "free-form" and long lines)
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN_ONCE([UD_FC_FREEFORM],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK([for Fortran flag needed to accept free-form source],
	       [ac_cv_fc_freeform],
[ac_cv_fc_freeform=unknown
ac_fc_freeform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffree-form -FR -free -qfree=f90 -qfree -Mfree -Mfreeform \
	       -freeform "-f free" -8 +source=free -nfix --nfix -Free
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_freeform_FCFLAGS_save $ac_flag"
dnl Use @&t@ below to ensure that editors don't turn 8+ spaces into tab.
  AC_COMPILE_IFELSE([[
  program freeform
       ! FIXME: how to best confuse non-freeform compilers?
       print *, 'Hello ', &
     @&t@     'world.'
       end]],
		    [ac_cv_fc_freeform=$ac_flag; break])
done
${RM} -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_freeform_FCFLAGS_save
])
if test "x$ac_cv_fc_freeform" = xunknown; then
  m4_default([$2],
	     [AC_MSG_WARN([Fortran $FC does not accept free-form source], 77)])
else
  dnl Do not append to FCFLAGS
  dnl if test "x$ac_cv_fc_freeform" != xnone; then
  dnl   FCFLAGS="$FCFLAGS $ac_cv_fc_freeform"
  dnl fi
  if test "x$ac_cv_fc_freeform" = xnone; then
     ac_cv_fc_freeform=
  fi
  $1
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_FREEFORM

dnl
dnl steal from autoconf 2.69
dnl customized AC_FC_FIXEDFORM: to just get ac_cv_fc_fixedform without
dnl appending ac_cv_fc_fixedform to FCFLAGS
dnl
# AC_FC_FIXEDFORM([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept
# fixed-format source code, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using new extension) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
#
# The known flags are:
#       -ffixed-form: GNU g77, gfortran, g95
#             -fixed: Intel compiler (ifort), Sun compiler (f95)
#            -qfixed: IBM compiler (xlf*)
#            -Mfixed: Portland Group compiler
#         -fixedform: SGI compiler
#           -f fixed: Absoft Fortran
#      +source=fixed: HP Fortran
#    (-)-fix, -Fixed: Lahey/Fujitsu Fortran
#             -fixed: NAGWare
# Since compilers may accept fixed form based on file name extension,
# but users may want to use it with others as well, call AC_FC_SRCEXT
# with the respective source extension before calling this macro.
AC_DEFUN_ONCE([UD_FC_FIXEDFORM],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK([for Fortran flag needed to accept fixed-form source],
               [ac_cv_fc_fixedform],
[ac_cv_fc_fixedform=unknown
ac_fc_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffixed-form -fixed -qfixed -Mfixed -fixedform "-f fixed" \
               +source=fixed -fix --fix -Fixed
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_fixedform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[
C     This comment should confuse free-form compilers.
      program main
      end]],
                    [ac_cv_fc_fixedform=$ac_flag; break])
done
${RM} -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_fixedform_FCFLAGS_save
])
if test "x$ac_cv_fc_fixedform" = xunknown; then
  m4_default([$2],
             [AC_MSG_WARN([Fortran does not accept fixed-form source], 77)])
  ac_cv_fc_fixedform=
else
  dnl Do not append to FCFLAGS
  dnl if test "x$ac_cv_fc_fixedform" != xnone; then
  dnl   FCFLAGS="$FCFLAGS $ac_cv_fc_fixedform"
  dnl fi
  if test "x$ac_cv_fc_fixedform" = xnone; then
     ac_cv_fc_fixedform=
  fi
  $1
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_FIXEDFORM

dnl Check if Fortran 77 compiler allows _8 modifier (a Fortran 90 feature)
dnl for integer*8 parameter
dnl NAG nagfor requires a modifier but does not like _8
dnl xlf is OK with _8 but when none is used it strangely passes the
dnl compilation with a warning message
dnl
AC_DEFUN([UD_FC_CONSTANT_MODIFIER],[
    AC_CACHE_CHECK([Fortran compiler treating constant modifier], [ac_cv_fc_constant_modifier],
    [AC_LANG_PUSH([Fortran 77])
        AC_COMPILE_IFELSE([[
         program main
         integer*8  nf_fill_uint
         integer*8  nf_fill_int64
         parameter (nf_fill_uint  = 4294967295_8)
         parameter (nf_fill_int64 = -9223372036854775806_8)
         end]],
        [ac_cv_fc_constant_modifier=8],
        [AC_COMPILE_IFELSE([[
         program main
         integer*8  nf_fill_uint
         integer*8  nf_fill_int64
         parameter (nf_fill_uint  = 4294967295)
         parameter (nf_fill_int64 = -9223372036854775806)
         end]],
        [ac_cv_fc_constant_modifier=none],
        [AC_COMPILE_IFELSE([[
         program main
         integer, parameter :: EightByteInt = selected_int_kind(18)
         integer*8  nf_fill_uint
         integer*8  nf_fill_int64
         parameter (nf_fill_uint  = 4294967295_EightByteInt)
         parameter (nf_fill_int64 = -9223372036854775806_EightByteInt)
         end]],
        [ac_cv_fc_constant_modifier=EightByteInt],
        [AC_MSG_ERROR([no appropriate modifier found])])
        ])
        ])
    ])
    AC_LANG_POP([Fortran 77])
])

dnl
dnl Borrowed macros from MPICH: aclocal_f77.m4
dnl

dnl
dnl Check to see if a C program can be linked when using the libraries
dnl needed by C programs
dnl
AC_DEFUN([PAC_PROG_FC_CHECK_FCLIBS],[
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_MSG_CHECKING([whether $CC links with FCLIBS found by autoconf])
AC_LANG_PUSH([C])
# Create a simple C program for the tests.
AC_LANG_CONFTEST([
    AC_LANG_PROGRAM([],[int a;])
])
# Try to link a C program with all of these libraries
saved_LIBS="$LIBS"
LIBS="$FCLIBS $saved_LIBS"
AC_LINK_IFELSE([],[
    AC_MSG_RESULT([yes])
],[
    AC_MSG_RESULT([no])
    AC_MSG_CHECKING([for which libraries can be used])
    pac_ldirs=""
    pac_libs=""
    pac_other=""
    for name in $FCLIBS ; do
        case $name in
        -l*) pac_libs="$pac_libs $name"   ;;
        -L*) pac_ldirs="$pac_ldirs $name" ;;
          *) pac_other="$pac_other $name" ;;
        esac
    done
    keep_libs=""
    for name in $pac_libs ; do
        LIBS="$saved_LIBS $pac_ldirs $pac_other $name"
        AC_LINK_IFELSE([],[
            keep_libs="$keep_libs $name"
        ])
    done
    AC_MSG_RESULT($keep_libs)
    FCLIBS="$pac_ldirs $pac_other $keep_libs"
])
LIBS="$saved_LIBS"
rm -f conftest.$ac_ext
AC_LANG_PUSH([C])
])

AC_DEFUN([UD_CHECK_F77_GNU_INT],
[
    AC_MSG_CHECKING([for Fortran 77 GNU intrinsic INT])
    AC_LANG_PUSH([Fortran 77])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
           program main
           implicit none
           real r
           integer*1 i1
           integer*2 i2
           integer   i4
           integer*8 i8
           r = 12.34
           i1 = INT(r, 1)
           i2 = INT(r, 2)
           i4 = INT(r)
           i8 = INT(r, 8)
           end
       ])],
       [ac_cv_f77_gnu_int="yes"],
       [ac_cv_f77_gnu_int="no"]
    )
    AC_LANG_POP([Fortran 77])
    AC_MSG_RESULT([$ac_cv_f77_gnu_int])
])

AC_DEFUN([UD_CHECK_F77_INT1],
[
    AC_MSG_CHECKING([for Fortran 77 intrinsic INT1])
    AC_LANG_PUSH([Fortran 77])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
           program main
           implicit none
           real r
           integer*1 i1
           r = 12.34
           i1 = INT1(r)
           end
       ])],
       [ac_cv_f77_int1="yes"],
       [ac_cv_f77_int1="no"]
    )
    AC_LANG_POP([Fortran 77])
    AC_MSG_RESULT([$ac_cv_f77_int1])
])

AC_DEFUN([UD_CHECK_F77_INT2],
[
    AC_MSG_CHECKING([for Fortran 77 intrinsic INT2])
    AC_LANG_PUSH([Fortran 77])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
           program main
           implicit none
           real r
           integer*2 i2
           r = 12.34
           i2 = INT2(r)
           end
       ])],
       [ac_cv_f77_int2="yes"],
       [ac_cv_f77_int2="no"]
    )
    AC_LANG_POP([Fortran 77])
    AC_MSG_RESULT([$ac_cv_f77_int2])
])

AC_DEFUN([UD_CHECK_F77_INT8],
[
    AC_MSG_CHECKING([for Fortran 77 intrinsic INT8])
    AC_LANG_PUSH([Fortran 77])
    AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([
           program main
           implicit none
           real r
           integer*8 i8
           r = 12.34
           i8 = INT8(r)
           end
       ])],
       [ac_cv_f77_int8="yes"],
       [ac_cv_f77_int8="no"]
    )
    AC_LANG_POP([Fortran 77])
    AC_MSG_RESULT([$ac_cv_f77_int8])
])

# _ACX_FC_MISMATCH([ACTION-IF-SUCCESS],
#                  [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to allow routines to be called with different
# argument types. The result is either "unknown", or the actual compiler flag
# required to downgrade consistency checking of procedure argument lists, which
# may be an empty string.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_mismatch_flag variable.
#
# Known flags:
# NAGWare: -mismatch
#
AC_DEFUN([_ACX_FC_MISMATCH],
  [_AC_FORTRAN_ASSERT()dnl
   m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_mismatch_flag])dnl
   AC_MSG_CHECKING([for _AC_LANG compiler flag needed to allow routines to dnl
be called with different argument types])
   AC_CACHE_VAL([acx_cache_var],
     [acx_cache_var=unknown
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      AC_LANG_CONFTEST([AC_LANG_PROGRAM([], [[      implicit none
      integer a
      real b
      character c
      call foo1(a)
      call foo1(b)
      call foo1(c)]])])
      for acx_flag in '' -mismatch; do
        _AC_LANG_PREFIX[]FLAGS="${acx_save_[]_AC_LANG_PREFIX[]FLAGS} $acx_flag"
        AC_COMPILE_IFELSE([], [acx_cache_var=$acx_flag])
        test "x$acx_cache_var" != xunknown && break
      done
      rm -f conftest.$ac_ext
      _AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
allow routines to be called with different argument types])])], [$1])
   m4_popdef([acx_cache_var])])

# ACX_FC_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
AC_DEFUN([ACX_FC_MISMATCH],
  [AC_LANG_PUSH([Fortran])_ACX_FC_MISMATCH($@)AC_LANG_POP([Fortran])])

# ACX_F77_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
AC_DEFUN([ACX_F77_MISMATCH],
  [AC_LANG_PUSH([Fortran 77])_ACX_FC_MISMATCH($@)AC_LANG_POP([Fortran 77])])

# _AC_PROG_FC_V
# -------------
# Determine the flag that causes the Fortran compiler to print
# information of library and object files (normally -v)
# Needed for _AC_FC_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: (-)-verbose, xlf: -V, Fujitsu: -###)
# Fujitsu accepts --verbose and passes it to the linker, which doesn't yield
# the desired result. Therefore test for -### before testing for --verbose.
# -------------
# This macro is used by AC_F77_LIBRARY_LDFLAGS and AC_FC_LIBRARY_LDFLAGS.
# We need to overload it to allow for additional possible results:
# NAG: -Wl,-v
AC_DEFUN([_AC_PROG_FC_V],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([how to get verbose linking output from $[]_AC_FC[]],
                [ac_cv_prog_[]_AC_LANG_ABBREV[]_v],
[AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_[]_AC_LANG_ABBREV[]_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose -V -\#\#\# --verbose -Wl,-v; do
  _AC_PROG_FC_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_[]_AC_LANG_ABBREV[]_v_output; do
     case $ac_arg in
	[[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
	  ac_cv_prog_[]_AC_LANG_ABBREV[]_v=$ac_verb
	  break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_[]_AC_LANG_ABBREV[]_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $[]_AC_FC[]])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_FC_V

dnl Check if Fortran 77 compiler is pgf77
dnl According to pgf77 manual the command-line option to should version is -V
dnl
dnl % pgf77 -V
dnl
dnl pgf77 16.9-0 64-bit target on x86-64 Linux -tp p7
dnl The Portland Group - PGI Compilers and Tools
dnl Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
dnl
dnl % pgfortran -V
dnl
dnl pgfortran 16.9-0 64-bit target on x86-64 Linux -tp p7
dnl The Portland Group - PGI Compilers and Tools
dnl Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
dnl
dnl Note the checking below may be obsolete
dnl
AC_DEFUN([UD_CHECK_F77_IS_PGF77],[
   if test "x$1" = x ; then
      TEST_F77=$F77
   else
      TEST_F77=$1
   fi
   AC_MSG_CHECKING([if Fortran 77 compiler is PGI pgf77])
   ac_cv_f77_is_PGF77=no
   eval $TEST_F77 -V </dev/null >& conftest.ver
   ac_F77_VENDOR=`head -c 5 conftest.ver`
   if test "x${ac_F77_VENDOR}" = xPGF77 ; then
      ac_cv_f77_is_PGF77=yes
   fi
   ${RM} -f conftest.ver
   unset ac_F77_VENDOR
   unset TEST_F77
   AC_MSG_RESULT([$ac_cv_f77_is_PGF77])
])
