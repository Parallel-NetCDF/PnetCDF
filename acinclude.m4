dnl $Id$
dnl UD macros for PnetCDF configure


dnl Convert a string to all uppercase.
dnl
define([uppercase],
[translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl
dnl Check for an m4(1) preprocessor utility.
dnl
AC_DEFUN([UD_PROG_M4],
[
    dnl AS_MESSAGE([checking for m4 preprocessor...])
    case "${M4-unset}" in
	unset) AC_CHECK_PROGS(M4, m4 gm4, []) ;;
	*) AC_CHECK_PROGS(M4, $M4 m4 gm4, []) ;;
    esac
    if test -z "$M4" ; then
       AC_MSG_ERROR("m4 utility program is required by PnetCDF")
    fi
    AC_MSG_CHECKING(m4 additional flags)
    case "${M4FLAGS-unset}" in
	unset) dnl test if M4 runs fine without option -B10000
               `${M4} /dev/null > conftest.err 2>&1`
               ac_cv_m4_stdout=`cat conftest.err`
               if test "x$ac_cv_m4_stdout" != x; then
                  M4FLAGS=-B10000
               fi
               ${RM} -f conftest.err
               ;;
    esac
    if test "x$M4FLAGS" = x; then
       AC_MSG_RESULT(none needed)
    else
       AC_MSG_RESULT($M4FLAGS)
    fi
    M4FFLAGS=`echo $M4FLAGS | sed 's/-s *//g'`
    dnl AC_MSG_NOTICE(M4FFLAGS=$M4FFLAGS)
    AC_SUBST(M4FLAGS)
    AC_SUBST(M4FFLAGS)
])

dnl
dnl Check for an ar(1) utility.
dnl
AC_DEFUN([UD_PROG_AR],
[
    dnl AS_MESSAGE([checking for ar utility...])
    case "${AR-unset}" in
	unset) AC_CHECK_PROGS(AR, ar, ar) ;;
	*) AC_CHECK_PROGS(AR, $AR ar, ar) ;;
    esac
    AC_MSG_CHECKING(ar flags)
    case "${ARFLAGS-unset}" in
	unset) ARFLAGS=cru ;;
    esac
    AC_MSG_RESULT($ARFLAGS)
    AC_SUBST(ARFLAGS)
])

dnl
dnl Check for an nm(1) utility.
dnl
AC_DEFUN([UD_PROG_NM],
[
    dnl AS_MESSAGE([checking for nm utility...])
    case "${NM-unset}" in
	unset) AC_CHECK_PROGS(NM, nm, nm) ;;
	*) AC_CHECK_PROGS(NM, $NM nm, nm) ;;
    esac
    AC_MSG_CHECKING(nm flags)
    case "${NMFLAGS-unset}" in
	unset) NMFLAGS= ;;
    esac
    AC_MSG_RESULT($NMFLAGS)
    AC_SUBST(NMFLAGS)
])

dnl
dnl Set the top-level source-directory.
dnl
AC_DEFUN([UD_SRCDIR],
[
    AC_MSG_CHECKING(for top-level source-directory)
    SRCDIR=`(cd $srcdir && pwd)`
    AC_MSG_RESULT($SRCDIR)
    AC_SUBST(SRCDIR)
])

dnl
dnl Check for a Standard C compiler.  Prefer a native one over the
dnl GNU one to reduce the chance that the environment variable LIBS
dnl will have to be set to reference the GNU C runtime library.
dnl
AC_DEFUN([UD_PROG_CC],
[
    # Because we must have a C compiler, we treat an unset CC
    # the same as an empty CC.
    case "${CC}" in
	'')
	    case `uname` in
		ULTRIX)
		    # The native ULTRIX C compiler isn't standard.
		    ccs='gcc cc'
		    ;;
		*)
		    # xlc is before c89 because AIX's sizeof(long long)
		    # differs between the two.
		    #
		    ccs='xlc c89 acc cc gcc'
		    ;;
	    esac
	    for cc in $ccs; do
		AC_CHECK_PROG(CC, $cc, $cc)
		case "$CC" in
		    '') ;;
		    *)  break
			;;
		esac
	    done
	    case "${CC}" in
		'')
		    AC_MSG_ERROR("Could not find C compiler")
		    ;;
	    esac
	    ;;
    esac
    #
    # On some systems, a discovered compiler nevertheless won't
    # work (due to licensing, for example); thus, we check the
    # compiler with a test program.
    #
    AC_MSG_CHECKING([C compiler "$CC"])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],
                      [AC_MSG_RESULT(works)],
                      [AC_MSG_RESULT(failed to compile test program)])
    AC_SUBST(CC)
    case "$CC" in
	*gcc*)
	    GCC=yes		# Expected by autoconf(1) macros
	    ;;
    esac
    case `uname -sr` in
	'HP-UX A.09'*)
	    AC_DEFINE(_HPUX_SOURCE)
	    ;;
    esac
])

dnl
dnl Check for a C++ compiler.  Prefer a native one over the
dnl GNU one to reduce the chance that the environment variable LIBS
dnl will have to be set to reference the GNU C runtime library.
dnl
AC_DEFUN([UD_PROG_CXX],
[
    case "${CXX-unset}" in
	unset)
	    case `uname` in
		AIX)
		    preferred_cxx='xlC'
		    ;;
	    esac
	    possible_cxxs="${preferred_cxx} CC cxx c++ g++ gcc"
	    ;;
	'') AC_MSG_WARN("Empty CXX variable")
	    possible_cxxs=
	    ;;
	*)  possible_cxxs=$CXX
	    ;;
    esac
    case "${possible_cxxs}" in
	'') CXX=
	    ;;
	*)  AC_LANG_PUSH([C++])
	    for cxx in $possible_cxxs; do
		AC_CHECK_PROG(CXX, $cxx, $cxx)
		case "$CXX" in
		    '') ;;
		    *)  # On some systems, a discovered compiler nevertheless
			# won't work (because it's a script to a non-existent
			# executable, for example); thus, we check the
			# compiler with a test program.  We also test
			# for <iostream.h> and the standard C++ library
			# because we need these to work.
			#
			AC_MSG_CHECKING(C++ compiler \"$CXX\")
			AC_RUN_IFELSE([AC_LANG_SOURCE([[
				#include <iostream.h>
				int main() {
				    cout << "";
				    return 0;
				}
			    ]])],[
				AC_MSG_RESULT(works)
				break
			    ],[
				AC_MSG_WARN($CXX failed on test program)
				CXX=
				unset ac_cv_prog_CXX
			    ],[])
			;;
		esac
	    done
	    AC_LANG_POP([C++])
	    case "${CXX}" in
		'') AC_MSG_WARN("Could not find working C++ compiler")
		    AC_MSG_WARN(Setting CXX to the empty string)
		    ;;
	    esac
	    ;;
    esac
    case "${CXX}" in
	'') AC_MSG_WARN(The C++ interface will not be built)
	    ;;
    esac
    AC_SUBST(CXX)
    case `uname` in
	'HP-UX A.09'*)
	    AC_DEFINE(_HPUX_SOURCE)
	    ;;
    esac
])


dnl
dnl like AC_TYPE_LONG_DOUBLE, except checks for 'long long'
dnl
AC_DEFUN([UD_C_LONG_LONG],
[AC_MSG_CHECKING(for long long)
AC_CACHE_VAL(ac_cv_c_long_long,
[if test "$GCC" = yes; then
  ac_cv_c_long_long=yes
else
AC_RUN_IFELSE([AC_LANG_SOURCE([[int main() {
long long foo = 0;
return(sizeof(long long) < sizeof(long)); }]])],
[ac_cv_c_long_long=yes],
[ac_cv_c_long_long=no],[:])
fi])dnl
AC_MSG_RESULT($ac_cv_c_long_long)
if test $ac_cv_c_long_long = yes; then
  AC_DEFINE(HAVE_LONG_LONG)
fi
])

dnl
dnl UD_CHECK_IEEE
dnl If the 'double' is not an IEEE double
dnl or the 'float' is not and IEEE single,
dnl define NO_IEEE_FLOAT
dnl
AC_DEFUN([UD_CHECK_IEEE],
[
AC_MSG_CHECKING(for IEEE floating point format)
AC_RUN_IFELSE([AC_LANG_SOURCE([[#ifndef NO_FLOAT_H
#include <float.h>
#endif

#define EXIT_NOTIEEE	1
#define EXIT_MAYBEIEEE	0

int
main()
{
#if	defined(FLT_RADIX)	&& FLT_RADIX != 2
		return EXIT_NOTIEEE;
#elif	defined(DBL_MAX_EXP)	&& DBL_MAX_EXP != 1024
		return EXIT_NOTIEEE;
#elif	defined(DBL_MANT_DIG)	&& DBL_MANT_DIG != 53
		return EXIT_NOTIEEE;
#elif 	defined(FLT_MAX_EXP)	&& !(FLT_MAX_EXP == 1024 || FLT_MAX_EXP == 128)
		return EXIT_NOTIEEE;
#elif	defined(FLT_MANT_DIG)	&& !(FLT_MANT_DIG == 53 || FLT_MANT_DIG == 24)
		return EXIT_NOTIEEE;
#else
	/* (assuming eight bit char) */
	if(sizeof(double) != 8)
		return EXIT_NOTIEEE;
	if(!(sizeof(float) == 4 || sizeof(float) == 8))
		return EXIT_NOTIEEE;

	return EXIT_MAYBEIEEE;
#endif
}]])],[ac_cv_c_ieeefloat=yes],[ac_cv_c_ieeefloat=no],[:])
AC_MSG_RESULT($ac_cv_c_ieeefloat)
if test x$ac_cv_c_ieeefloat = xno; then
  AC_DEFINE(NO_IEEE_FLOAT)
fi
])

dnl Check for utility for generating makefile dependencies.
dnl Should only be used at the UPC.
dnl
AC_DEFUN([UD_PROG_CC_MAKEDEPEND],
[
    AC_MSG_CHECKING(how to make dependencies)
    case `uname -s` in
	IRIX*|OSF1)
	    CC_MAKEDEPEND='cc -M'
	    ;;
	SunOS)
	    case `uname -r` in
		4*)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		5*|*)
		    CC_MAKEDEPEND='cc -xM'
		    ;;
	    esac
	    ;;
	ULTRIX)
	    case `uname -m` in
		RISC)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		VAX)	# Can't handle prototypes in netcdf.h
		    ;;
	    esac
	    ;;
	AIX)	# Writes to .u files rather than standard out
	    ;;
	HP-UX)	# Writes escaped newlines to standard error
	    ;;
    esac
    case "${CC_MAKEDEPEND}" in
	'')
	    CC_MAKEDEPEND=false
	    ;;
    esac
    AC_MSG_RESULT($CC_MAKEDEPEND)
    AC_SUBST(CC_MAKEDEPEND)
])


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


dnl Setup for making a manual-page database.
dnl
AC_DEFUN([UD_MAKEWHATIS],
[
    #
    # NB: We always want to define WHATIS to prevent the
    # $(MANDIR)/$(WHATIS) make(1) target from being just $(MANDIR)/ and
    # conflicting with the (directory creation) target with the same name.
    #
    WHATIS=whatis
    case `uname -sr` in
	BSD/OS*|FreeBSD*)
	    # Can't generate a user-database -- only /usr/share/man/whatis.db.
	    MAKEWHATIS_CMD=
	    ;;
	'IRIX64 6.5'|'IRIX 6.5')
	    MAKEWHATIS_CMD='/usr/lib/makewhatis -M $(MANDIR) $(MANDIR)/whatis'
	    ;;
	'IRIX 6'*)
	    # Can't generate a user-database.
	    MAKEWHATIS_CMD=
	    ;;
	HP-UX*)
	    # Can't generate a user-database -- only /usr/lib/whatis.
	    MAKEWHATIS_CMD=
	    ;;
	'Linux '*)
	    # /usr/sbin/makewhatis doesn't work
	    MAKEWHATIS_CMD=
	    ;;
	ULTRIX*)
	    # Can't generate a user-database -- only /usr/lib/whatis.
	    MAKEWHATIS_CMD=
	    ;;
	*)
	    if test -r /usr/man/windex; then
		WHATIS=windex
	    fi
	    AC_CHECK_PROGS(prog, catman makewhatis /usr/lib/makewhatis)
	    case "$prog" in
		*catman*)
		    MAKEWHATIS_CMD=$prog' -w -M $(MANDIR)'
		    ;;
		*makewhatis*)
		    MAKEWHATIS_CMD=$prog' $(MANDIR)'
		    ;;
	    esac
	    ;;
    esac
    AC_SUBST(WHATIS)
    AC_SUBST(MAKEWHATIS_CMD)
    AC_MSG_CHECKING(for manual-page index command)
    AC_MSG_RESULT($MAKEWHATIS_CMD)
])


dnl Check for the math library.
dnl
AC_DEFUN([UD_CHECK_LIB_MATH],
[
    dnl AS_MESSAGE([checking for math library...])
    case "${MATHLIB}" in
	'')
	    AC_CHECK_LIB(c, floor, MATHLIB=,
		    AC_CHECK_LIB(m, floor, MATHLIB=-lm, MATHLIB=))
	    ;;
	*)
	    AC_MSG_RESULT($MATHLIB (user defined))
	    ;;
    esac
    AC_SUBST(MATHLIB)
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

# AX_C_FLOAT_WORDS_BIGENDIAN# added by:
#   Warren Turkal <wt@penguintechs.org>
#
# Copyright © 2006 Daniel Amelang <dan@amelang.net>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and
# this notice are preserved.
#
# This macro will detect if double variables are words packed in big endian
# order while the bits in the words are arranged in little endian order. This
# macro was added to support the ARM architecture. The FLOAT_WORDS_BIGENDIAN
# macro will be set to 1 if the word order is big endian. If the word order is
# not big endian, FLOAT_WORDS_BIGENDIAN will be not be set.
AC_DEFUN([AX_C_FLOAT_WORDS_BIGENDIAN],
  [AC_CACHE_CHECK(whether float word ordering is bigendian,
                  ax_cv_c_float_words_bigendian, [

ax_cv_c_float_words_bigendian=unknown
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[

double d = 9090423496703681033747047890550501147621169273561563201479712084405348
886581669527372346909785805625751702019124748742951693213050356065000232756451757
0778480236724525140520121371739201496540132640109977779420565776568942592.0;

]])], [

if grep noonsees conftest.$ac_objext >/dev/null ; then
  ax_cv_c_float_words_bigendian=yes
fi
if grep seesnoon conftest.$ac_objext >/dev/null ; then
  if test "$ax_cv_c_float_words_bigendian" = unknown; then
    ax_cv_c_float_words_bigendian=no
  else
    ax_cv_c_float_words_bigendian=unknown
  fi
fi

])])

case $ax_cv_c_float_words_bigendian in
  yes)
    m4_default([$1],
      [AC_DEFINE([FLOAT_WORDS_BIGENDIAN], 1,
                 [Define to 1 if your system stores words within floats
                  with the most significant word first])]) ;;
  no)
    $2 ;;
  *)
    m4_default([$3],
      [AC_MSG_ERROR([

Unknown float word ordering. You need to manually preset
ax_cv_c_float_words_bigendian=no (or yes) according to your system.

    ])]) ;;
esac

])# AX_C_FLOAT_WORDS_BIGENDIAN

dnl Find the full path of a header file
dnl
dnl UD_CHECK_HEADER_PATH(file, [action-if-found], [action-if-not-found])
dnl Example:
dnl UD_CHECK_HEADER_PATH([math.h])
dnl AC_MSG_NOTICE([ac_cv_header_path_math_h=$ac_cv_header_path_math_h])
dnl
dnl
AC_DEFUN([UD_CHECK_HEADER_PATH],
[
    AS_VAR_PUSHDEF([ac_Path], [ac_cv_header_path_$1])dnl
    AC_CACHE_CHECK(
       [for full path of header file \"$1\"], [ac_Path],
       [AC_PREPROC_IFELSE(
           [AC_LANG_PROGRAM([[#include <$1>]])],
           [AS_VAR_SET([ac_Path], [`sed -n '/\.h"/s/.*"\(.*\)".*/\1/p' conftest.i | grep -m 1 $1`])],
           [AC_MSG_RESULT([not found])]
       )])
    AS_VAR_SET_IF([ac_Path], [$2], [$3])
    dnl eval AS_TR_SH([ac_cv_header_path_$1])=$ac_Path
    AS_VAR_POPDEF([ac_Path])dnl
])

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

dnl Check if MPI compiler is GNU
dnl According to gcc manual the command-line option to show version is --version
dnl
dnl % gcc --version
dnl gcc (Ubuntu 4.8.4-2ubuntu1~14.04.4) 4.8.4
dnl Copyright (C) 2013 Free Software Foundation, Inc.
dnl This is free software; see the source for copying conditions.  There is NO
dnl warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_GCC],[
    AC_CACHE_CHECK([if MPI C compiler is GNU gcc], [ac_cv_mpicc_is_GCC],
    [ac_cv_mpicc_is_GCC=no
     ac_MPICC_VER=`$MPICC --version`
     ac_MPICC_VENDOR=`echo $ac_MPICC_VER | cut -s -d' ' -f1`
     if test "x${ac_MPICC_VENDOR}" = xgcc ; then
        ac_cv_mpicc_is_GCC=yes
     fi
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is CLANG
dnl According to clang manual the command-line option to show version is --version
dnl
dnl % clang --version
dnl Ubuntu clang version 3.4-1ubuntu3 (tags/RELEASE_34/final) (based on LLVM 3.4)
dnl Target: x86_64-pc-linux-gnu
dnl Thread model: posix
dnl
dnl or
dnl
dnl clang version 3.4.2 (tags/RELEASE_34/dot2-final)
dnl Target: x86_64-redhat-linux-gnu
dnl Thread model: posix
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_CLANG],[
    AC_CACHE_CHECK([if MPI C compiler is Clang], [ac_cv_mpicc_is_CLANG],
    [ac_cv_mpicc_is_CLANG=no
     ac_MPICC_VER=`$MPICC --version`
     ac_MPICC_VENDOR=`echo $ac_MPICC_VER | ${GREP} -w clang`
     if test "x${ac_MPICC_VENDOR}" != x ; then
        ac_cv_mpicc_is_CLANG=yes
     fi
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is Intel icc
dnl According to icc manual the command-line option to show version is --version
dnl
dnl % icc --version
dnl icc (ICC) 17.0.0 20160721
dnl Copyright (C) 1985-2016 Intel Corporation.  All rights reserved.
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_ICC],[
    AC_CACHE_CHECK([if MPI C compiler is Intel icc], [ac_cv_mpicc_is_ICC],
    [ac_cv_mpicc_is_ICC=no
     ac_MPICC_VER=`$MPICC --version`
     ac_MPICC_VENDOR=`echo $ac_MPICC_VER | cut -s -d' ' -f1`
     if test "x${ac_MPICC_VENDOR}" = xicc ; then
        ac_cv_mpicc_is_ICC=yes
     fi
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is Fujitsu fccpx based
dnl According to mpifccpx manual the command-line option to show version is
dnl -showme
dnl
dnl % mpifccpx --showme
dnl /opt/FJSVtclang/GM-1.2.0-24/bin/fccpx -Kident_mpi -mt ...
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_FCCPX],[
    AC_CACHE_CHECK([if MPI C compiler is Fujitsu fccpx], [ac_cv_mpicc_is_FCCPX],
    [ac_cv_mpicc_is_FCCPX=no
     ac_MPICC_VENDOR=`$MPICC --showme 2> /dev/null | cut -d' ' -f1 | xargs -r basename`
     if test "x${ac_MPICC_VENDOR}" = xfccpx ; then
        ac_cv_mpicc_is_FCCPX=yes
     fi
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is IBM XL based
dnl According to xlc manual the command-line option to show version is
dnl
dnl % xlc -qversion
dnl IBM XL C/C++ for Blue Gene, V12.1
dnl Version: 12.01.0000.0011
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_XLC],[
    AC_CACHE_CHECK([if MPI C compiler is IBM XLC], [ac_cv_mpicc_is_XLC],
    [ac_cv_mpicc_is_XLC=no
     ac_MPICC_VER=`$MPICC -qversion >& conftest.ver`
     ac_MPICC_VENDOR=`head -c 6 conftest.ver`
     if test "x${ac_MPICC_VENDOR}" = "xIBM XL" ; then
        ac_cv_mpicc_is_XLC=yes
     fi
     ${RM} -f conftest.ver
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is pgcc based
dnl According to pgcc manual the command-line option to show version is -V
dnl
dnl % pgcc -V
dnl
dnl pgcc 16.9-0 64-bit target on x86-64 Linux -tp p7
dnl The Portland Group - PGI Compilers and Tools
dnl Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_PGCC],[
    AC_CACHE_CHECK([if MPI C compiler is PGI pgcc], [ac_cv_mpicc_is_PGCC],
    [ac_cv_mpicc_is_PGCC=no
     ac_MPICC_VER=`$MPICC -V -c 2> /dev/null`
     ac_MPICC_VENDOR=`echo $ac_MPICC_VER | cut -s -d' ' -f1`
     if test "x${ac_MPICC_VENDOR}" = xpgcc ; then
        ac_cv_mpicc_is_PGCC=yes
     fi
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check if MPI compiler is Oracle Solaris Studio
dnl According to cc manual the command-line option to show version is -V
dnl
dnl % cc -V
dnl cc: Sun C 5.13 Linux_i386 2014/10/20
dnl
AC_DEFUN([UD_CHECK_MPICC_IS_SOLARIS],[
    AC_CACHE_CHECK([if MPI C compiler is Solaris cc], [ac_cv_mpicc_is_SOLARIS],
    [ac_cv_mpicc_is_SOLARIS=no
     ac_MPICC_VER="$($MPICC -V 2>&1)"
     ac_MPICC_VENDOR=`echo $ac_MPICC_VER | cut -s -d' ' -f2`
     if test "x${ac_MPICC_VENDOR}" = xSun ; then
        ac_cv_mpicc_is_SOLARIS=yes
     fi
     unset ac_MPICC_VER
     unset ac_MPICC_VENDOR
    ])
])

dnl Check MPI C compiler base
dnl
AC_DEFUN([UD_CHECK_MPICC_BASE_VENDOR],[
   AC_MSG_CHECKING([MPI C compiler base])
   ac_cv_mpicc_base_vendor=
   # Check GCC
   ac_MPICC_VER="$($MPICC --version 2>&1)"
   ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w gcc`
   # AC_MSG_NOTICE(GCC ac_MPICC_VER=$ac_MPICC_VER)
   if test "x${ac_MPICC_VER}" != x ; then
      ac_cv_mpicc_base_vendor="GCC"
   else
      # Check CLANG
      ac_MPICC_VER="$($MPICC --version 2>&1)"
      ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w clang`
      # AC_MSG_NOTICE(clang ac_MPICC_VER=$ac_MPICC_VER)
      if test "x${ac_MPICC_VER}" != x ; then
         ac_cv_mpicc_base_vendor="CLANG"
      else
         # Check Intel C
         ac_MPICC_VER="$($MPICC --version 2>&1)"
         # grep keyword Intel instead of "icc"
         ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w Intel`
         # AC_MSG_NOTICE(icc ac_MPICC_VER=$ac_MPICC_VER)
         if test "x${ac_MPICC_VER}" != x ; then
            ac_cv_mpicc_base_vendor="ICC"
         else
            # Check XLC
            ac_MPICC_VER="$($MPICC -qversion 2>&1)"
            ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} "IBM XL C"`
            # AC_MSG_NOTICE(XLC ac_MPICC_VER=$ac_MPICC_VER)
            if test "x${ac_MPICC_VER}" != x ; then
               ac_cv_mpicc_base_vendor="XLC"
            else
               # Check PGCC
               ac_MPICC_VER="$($MPICC -V -c 2>&1)"
               # grep keyword PGI instead of "pgcc"
               ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w PGI`
               # AC_MSG_NOTICE(pgcc ac_MPICC_VER=$ac_MPICC_VER)
               if test "x${ac_MPICC_VER}" != x ; then
                  ac_cv_mpicc_base_vendor="PGCC"
               else
                  # Check SOLARIS
                  ac_MPICC_VER="$($MPICC -V 2>&1)"
                  ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w Sun`
                  # AC_MSG_NOTICE(Sun ac_MPICC_VER=$ac_MPICC_VER)
                  if test "x${ac_MPICC_VER}" != x ; then
                     ac_cv_mpicc_base_vendor="SOLARIS"
                  else
                     # Check FCCPX
                     ac_MPICC_VER="$($MPICC --showme 2>&1)"
                     ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w fccpx`
                     # AC_MSG_NOTICE(fccpx ac_MPICC_VER=$ac_MPICC_VER)
                     if test "x${ac_MPICC_VER}" != x ; then
                        ac_cv_mpicc_base_vendor="FCCPX"
                     else
                        # If just cc, check if it is a wrapper of GCC
                        ac_MPICC_VER="$($MPICC -v 2>&1)"
                        UD_MSG_DEBUG(GCC ac_MPICC_VER=$ac_MPICC_VER)
                        ac_MPICC_VER=`echo $ac_MPICC_VER | ${GREP} -w gcc`
                        if test "x${ac_MPICC_VER}" != x ; then
                           ac_cv_mpicc_base_vendor="GCC"
                        fi
                     fi
                  fi
               fi
            fi
         fi
      fi
   fi
   if test "x$ac_cv_mpicc_base_vendor" = x ; then
      AC_MSG_RESULT([unknown])
   else
      AC_MSG_RESULT([$ac_cv_mpicc_base_vendor])
   fi
])

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
AC_DEFUN([UD_CHECK_MPIF77_IS_PGF77],[
    AC_CACHE_CHECK([if MPI Fortran 77 compiler is PGI pgf77], [ac_cv_mpif77_is_PGF77],
    [ac_cv_mpif77_is_PGF77=no
     eval $MPIF77 -V </dev/null >& conftest.ver
     ac_MPIF77_VENDOR=`head -c 5 conftest.ver`
     if test "x${ac_MPIF77_VENDOR}" = xPGF77 ; then
        ac_cv_mpif77_is_PGF77=yes
     fi
     ${RM} -f conftest.ver
     unset ac_MPIF77_VENDOR
    ])
])

AC_DEFUN([UD_CXX_MACRO_FUNC],[
   AC_CACHE_CHECK([if C++ macro __func__ or __FUNCTION__ is defined], [ac_cv_cxx_macro_func],
   [ac_cv_cxx_macro_func=no
    ac_cv_cxx_macro_function=no
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <iostream>]],
                                       [[std::cout << __func__;]])],
                                       [ac_cv_cxx_macro_func=yes],[])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <iostream>]],
                                       [[std::cout << __FUNCTION__;]])],
                                       [ac_cv_cxx_macro_function=yes],[])
    AC_LANG_POP([C++])
   ])
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

dnl
dnl check the availability of one MPI executable in $2
dnl
dnl $2 can be a single command, This is the case when user set the environment.
dnl The variable may contain the executable name followed by zeor or more
dnl command-line options. In the latter case, we check the first string token,
dnl the command name, and ignore the rest command-line options. For example,
dnl UD_MPI_PATH_PROG([MPICC], [mpicc -O2])
dnl
dnl In addition, the first token of $2 may contain the full path of
dnl the command. For example, UD_MPI_PATH_PROG([MPICC], [/usr/bin/mpicc -O2])
dnl
AC_DEFUN([UD_MPI_PATH_PROG], [
   if test "x$2" = x ; then
      AC_MSG_ERROR("2nd argument cannot be NULL")
   fi

   dnl 1st token in $2 must be the program name, rests are command-line options
   ac_first_token=`echo $2 | cut -d" " -f1`
   ac_rest_tokens=`echo $2 | cut -d" " -s -f2-`
   UD_MSG_DEBUG(ac_first_token=$ac_first_token) dnl executable name
   UD_MSG_DEBUG(ac_rest_tokens=$ac_rest_tokens) dnl command-line option

   dnl First check if ac_first_token contain a full path
   dnl If yes, check, check if the file exists. Need not check MPI_INSTALL.
   ac_mpi_prog_path=`AS_DIRNAME(["$ac_first_token"])`
   if test "x$ac_mpi_prog_path" != "x." ; then
      AC_MSG_CHECKING([whether $ac_first_token exists and is executable])
      if test -x "$ac_first_token" ; then
         AC_MSG_RESULT([yes])
         $1="$2"
      else
         AC_MSG_RESULT([no])
         $1=
      fi
   else
      dnl ac_first_token does not contain a full path
      ac_mpi_prog_$1=
      if test "x$MPI_INSTALL" != x ; then
         dnl First, check if it can be found under $MPI_INSTALL, i.e.
         dnl --with-mpi is used on configure command line
         if test -d "${MPI_INSTALL}/bin" ; then
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}/bin])
            if test -x "$MPI_INSTALL/bin/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/bin/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         else
            dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}])
            if test -x "$MPI_INSTALL/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         fi
         if test "x$ac_mpi_prog_$1" != x ; then
            if test "x$ac_rest_tokens" != x ; then
               $1="${ac_mpi_prog_$1} $ac_rest_tokens"
            else
               $1=${ac_mpi_prog_$1}
            fi
         else
            $1=
         fi
      else
         dnl MPI_INSTALL is not set, i.e. --with-mpi is not used
         AC_PATH_PROG([ac_mpi_prog_$1], [$ac_first_token])
         if test "x$ac_mpi_prog_$1" != x ; then
            if test "x$ac_rest_tokens" != x ; then
               $1="${ac_mpi_prog_$1} $ac_rest_tokens"
            else
               $1=${ac_mpi_prog_$1}
            fi
         else
            $1=
         fi
      fi
   fi
])

dnl
dnl check the availability of a list of MPI executables
dnl
dnl Note $2 can be a list of executable commands to be searched with each
dnl command being the executable file name without command-line option. This is
dnl the case when user does not set the environment variable, for example
dnl MPICC, and we must search one from the candidate list. For example,
dnl UD_MPI_PATH_PROGS([MPICC], [mpicc mpixlc mpifccpx mpipgcc])
dnl
AC_DEFUN([UD_MPI_PATH_PROGS], [
   ac_mpi_prog_$1=
   if test "x$MPI_INSTALL" != x ; then
      UD_MSG_DEBUG(--with-mpi=$MPI_INSTALL is used)
      if test -d "${MPI_INSTALL}/bin" ; then
         UD_MSG_DEBUG(search $2 under $MPI_INSTALL/bin)
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL/bin])
      else
         dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
         UD_MSG_DEBUG(search $2 under $MPI_INSTALL)
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL])
      fi
   else
      UD_MSG_DEBUG(--with-mpi=$MPI_INSTALL is NOT used)
      UD_MSG_DEBUG(search $2 under $PATH)
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
   fi
   UD_MSG_DEBUG([ac_mpi_prog_$1=${ac_mpi_prog_$1}])
   if test "x${ac_mpi_prog_$1}" = x ; then
      dnl AC_CHECK_FILES fails when $2 is not found in cross compile
      dnl AC_CHECK_FILES([$2], [ac_mpi_prog_$1=$2])
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2], [], [/])
      dnl ac_first_token=`echo $2 | cut -d" " -f1`
      dnl UD_MSG_DEBUG(check first token $ac_first_token of $2)
      dnl if test -f $ac_first_token ; then
         dnl UD_MSG_DEBUG(use file $ac_first_token as it exits)
         dnl ac_mpi_prog_$1=$2
      dnl fi
   fi
   $1=${ac_mpi_prog_$1}
])

dnl Check for presence of an MPI constant.
dnl These could be enums, so we have to do compile checks.
AC_DEFUN([UD_HAS_MPI_CONST], [
   AC_MSG_CHECKING(whether MPI constant $1 is defined )
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <mpi.h>
          int dummy = $1;
      ])],
      [AC_MSG_RESULT(yes)
       AC_DEFINE(HAVE_$1, 1, available)
      ],
      [AC_MSG_RESULT(no)]
   )]
)

dnl Check for presence of an MPI datatype.
dnl These could be enums, so we have to do compile checks.
AC_DEFUN([UD_CHECK_MPI_DATATYPE], [
   AC_MSG_CHECKING(whether MPI datatype $1 is defined )
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <mpi.h>
          MPI_Datatype dummy = $1;
      ])],
      [ac_cv_CHECK_MPI_DATATYPE_$1=yes],
      [ac_cv_CHECK_MPI_DATATYPE_$1=no]
   )
   AC_MSG_RESULT($ac_cv_CHECK_MPI_DATATYPE_$1)
   if test "x$ac_cv_CHECK_MPI_DATATYPE_$1" = xyes; then
       AC_DEFINE(HAVE_DECL_$1, 1, available)
   fi
])

dnl Check if older Intel MPI C compiler (4.x) for issue of redefined SEEK_SET
dnl See https://software.intel.com/en-us/articles/intel-cluster-toolkit-for-linux-error-when-compiling-c-aps-using-intel-mpi-library-compilation-driver-mpiicpc
AC_DEFUN([UD_CHECK_MPI_CPP_SEEK_SET], [
   AC_MSG_CHECKING(whether MPI C++ compiler redefines SEEK_SET )
   CXX=${MPICXX}
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <stdio.h>
          #include <mpi.h>
          int main() { return 0; }
      ])],
      [ac_cv_CHECK_MPI_CPP_SEEK_SET=no],
      [ac_cv_CHECK_MPI_CPP_SEEK_SET=yes]
   )
   AC_MSG_RESULT([$ac_cv_CHECK_MPI_CPP_SEEK_SET])
   AC_LANG_POP(C++)
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

dnl
dnl Check how sed command handling in-place option -i and define SED_I
dnl
AC_DEFUN([UD_PROG_SED_I],
[
   AC_REQUIRE([AC_PROG_SED])
   AC_CACHE_CHECK([for sed handling option -i ], ac_cv_SED_I,[
   cat > conftest.sed_i <<EOF
   test str1
EOF
   ac_cv_err=`$SED -i '' -e 's|str1|str2|g' conftest.sed_i 2>&1`
   if test "x$ac_cv_err" = x ; then
      ac_cv_SED_I="$SED -i ''"
   else
      ac_cv_err=`sed -i'' -e 's|str1|str2|g' conftest.sed_i 2>&1`
      if test "x$ac_cv_err" = x ; then
         ac_cv_SED_I="$SED -i''"
      else
         AC_MSG_ERROR("No proper sed -i option found")
      fi
   fi
   AS_UNSET(ac_cv_err)])
   SED_I="$ac_cv_SED_I"
   AC_SUBST(SED_I)
   ${RM} -f conftest.sed_i
])

#
# To check whether MPI C compiler supports shared libraries, below three
# functions LT_AC_CHECK_SHLIB, LT_AH_CHECK_SHLIB, and LT_AC_LINK_SHLIB_IFELSE
# are from http://lists.gnu.org/archive/html/libtool/2004-10/msg00222.html
#
# Usage example:
#
#  LT_AC_CHECK_SHLIB(mpi, MPI_Init, [], [AC_MSG_ERROR([
#    -----------------------------------------------------------------------
#      Building shared libraries is requested, but the MPI C compiler:
#      "${MPICC}"
#      is not built with support of shared libraries. Abort.
#    -----------------------------------------------------------------------])])
#

# LT_AC_CHECK_SHLIB(LIBRARY, FUNCTION,
#                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#                   [OTHER-LIBRARIES])
# -----------------------------------------------------------
#
# Use a cache variable name containing both the library and function name,
# because the test really is for library $1 defining function $2, not
# just for library $1. Separate tests with the same $1 and different $2s
# may have different results.
#
# Note that using directly AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1_$2])
# is asking for troubles, since LT_AC_CHECK_SHLIB($lib, fun) would give
# lt_ac_cv_shlib_$lib_fun, which is definitely not what was meant. Hence
# the AS_LITERAL_IF indirection.
#
# FIXME: This macro is extremely suspicious. It DEFINEs unconditionally,
# whatever the FUNCTION, in addition to not being a *S macro.  Note
# that the cache does depend upon the function we are looking for.
#
AC_DEFUN([LT_AC_CHECK_SHLIB],
[m4_ifval([$3], , [LT_AH_CHECK_SHLIB([$1])])dnl
AS_LITERAL_IF([$1],
              [AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1_$2])],
              [AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1''_$2])])dnl
AC_CACHE_CHECK([for $2 in shared version of -l$1], ac_Lib,
[lt_ac_check_shlib_save_LIBS=$LIBS
LIBS="-l$1 $5 $LIBS"
LT_AC_LINK_SHLIB_IFELSE([AC_LANG_CALL([], [$2])],
                        [AS_VAR_SET(ac_Lib, yes)],
                        [AS_VAR_SET(ac_Lib, no)])
LIBS=$lt_ac_check_shlib_save_LIBS])
AS_IF([test AS_VAR_GET(ac_Lib) = yes],
      [m4_default([$3], [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_SHLIB$1))])],
      [$4])dnl
AS_VAR_POPDEF([ac_Lib])dnl
])# AC_CHECK_LIB

# LT_AH_CHECK_SHLIB(LIBNAME)
# ---------------------
m4_define([LT_AH_CHECK_SHLIB],
[AH_TEMPLATE(AS_TR_CPP(HAVE_SHLIB$1),
[Define to 1 if you have a shared version of the `]$1[' library (-l]$1[).])])


# LT_AC_LINK_SHLIB_IFELSE(LIBRARYCODE, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------
# Try to link LIBRARYCODE into a libtool library.
AC_DEFUN([LT_AC_LINK_SHLIB_IFELSE],
[m4_ifvaln([$1], [AC_LANG_CONFTEST([$1])])dnl
rm -rf $objdir
rm -f conftest.$ac_objext conftest.la
ac_ltcompile='./libtool --mode=compile $CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext -o conftest.lo >&AS_MESSAGE_LOG_FD'
ac_ltlink_la='./libtool --mode=link $CC -rpath `pwd` $CFLAGS $LDFLAGS -o libconftest.la conftest.lo $LIBS >&AS_MESSAGE_LOG_FD'
AS_IF([AC_TRY_EVAL([ac_ltcompile]) &&
       AC_TRY_EVAL([ac_ltlink_la]) &&
       AC_TRY_COMMAND([test -s libconftest.la])],
      [$2],
      [_AC_MSG_LOG_CONFTEST
      m4_ifvaln([$3], [$3])])
rm -rf $objdir
rm -f conftest* libconftest*[]dnl
])# _AC_LINK_IFELSE

# LT_MPI_CHECK_SHLIB
# -----------------------------------------------------------------
# Try to link an MPI program using libtool. This function is useful for
# detecting whether the MPI library is built with shared library support.
AC_DEFUN([LT_MPI_CHECK_SHLIB],[
   AC_MSG_CHECKING([whether MPI library is built with shared library support])
   AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <mpi.h>]], [[MPI_Init(0, 0);]])])
   $RM -rf $objdir conftest.$ac_objext conftest.la
   dnl RM must have -f option when calling libtool
   ac_RM_saved=${RM}
   if test "x$RM" = xrm || test "x$RM" = "x/bin/rm" ; then
      RM="$RM -f"
   fi
   ac_ltcompile='./libtool --mode=compile $MPICC -c $CFLAGS $CPPFLAGS conftest.$ac_ext -o conftest.lo >&AS_MESSAGE_LOG_FD'
   ac_ltlink_la='./libtool --mode=link $MPICC -rpath `pwd` $CFLAGS $LDFLAGS -o libconftest.la conftest.lo $LIBS >&AS_MESSAGE_LOG_FD'
   AS_IF([AC_TRY_EVAL([ac_ltcompile]) &&
       AC_TRY_EVAL([ac_ltlink_la]) &&
       AC_TRY_COMMAND([test -s libconftest.la])],
      [ac_cv_lt_mpi_check_shlib=yes],
      [ac_cv_lt_mpi_check_shlib=no])
   RM=$ac_RM_saved
   unset ac_RM_saved
   $RM -rf $objdir conftest* libconftest*
   AC_MSG_RESULT([$ac_cv_lt_mpi_check_shlib])
])# LT_MPI_CHECK_SHLIB

# CHECK_MPI_VERSION
# -----------------------------------------------------------------
# check MPI version and vendor
AC_DEFUN([CHECK_MPI_VERSION],[
   AC_REQUIRE([AX_COMPILER_VENDOR])
   AC_REQUIRE([AC_PROG_GREP])
   AC_MSG_CHECKING([MPI Standard version implemented])
   AC_COMPUTE_INT([mpi_version], [MPI_VERSION], [[#include <mpi.h>]])
   AC_COMPUTE_INT([mpi_subversion], [MPI_SUBVERSION], [[#include <mpi.h>]])
   if test "x$mpi_version" = x ; then
      AC_MSG_RESULT([information unavailable])
   else
      AC_MSG_RESULT([${mpi_version}.${mpi_subversion}])
   fi

   AC_CHECK_DECL([MPICH_VERSION],      [], [], [#include <mpi.h>])
   AC_CHECK_DECL([MPICH2_VERSION],     [], [], [#include <mpi.h>])
   AC_CHECK_DECL([OMPI_MAJOR_VERSION], [], [], [#include <mpi.h>])
   AC_CHECK_DECL([MVAPICH2_VERSION],   [], [], [#include <mpi.h>])
   AC_MSG_CHECKING([MPI vendor])

   # MACRO_FLAG is the CPP flag to show macro definitions. For GCC, Intel, and
   # PGI, it is -dM and prints to stdout. For Oracle Solaris Studio compiler,
   # it is -xdumpmacros=defs and prints to stderr
   MACRO_FLAG="-dM"
   if test "x$ax_cv_c_compiler_vendor" = xibm ; then
      MACRO_FLAG="-qshowmacros"
   fi

   saved_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $MACRO_FLAG"
   AC_PREPROC_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [])],
                     [`cp conftest.i saved_conftest.i`])
   CPPFLAGS=$saved_CPPFLAGS
   unset MACRO_FLAG

   # version_str=`${EGREP} -E 'MVAPICH2_VERSION|MPICH_VERSION|MPICH2_VERSION|OMPI_MAJOR_VERSION|OMPI_MINOR_VERSION|OMPI_RELEASE_VERSION' conftest.i`
   # echo a variable deletes any linefeeds in that variable, so we cannot use
   # `echo $version_str | ${GREP} MVAPICH2_VERSION`. Instead, we can use
   # version=`${GREP} MPICH_VERSION <<< "$version_str" | cut -d' ' -d'"' -f2`

   # Note MVAPICH2's mpi.h also defines MPICH_VERSION, so this check must be
   # done before MPICH.
   if test -f saved_conftest.i ; then
      if test "x$ac_cv_have_decl_MVAPICH2_VERSION" = xyes ; then
         mvapich2_version=`${GREP} MVAPICH2_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MVAPICH2 $mvapich2_version)
         unset mvapich2_version
      elif test "x$ac_cv_have_decl_MPICH_VERSION" = xyes ; then
         mpich_version=`${GREP} MPICH_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MPICH $mpich_version)
         unset mpich_version
      elif test "x$ac_cv_have_decl_MPICH2_VERSION" = xyes ; then
         mpich2_version=`${GREP} MPICH2_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MPICH2 $mpich2_version)
         unset mpich2_version
      elif test "x$ac_cv_have_decl_OMPI_MAJOR_VERSION" = xyes ; then
         # AC_COMPUTE_INT([OMPI_MAJOR], [OMPI_MAJOR_VERSION], [[#include <mpi.h>]])
         # AC_COMPUTE_INT([OMPI_MINOR], [OMPI_MINOR_VERSION], [[#include <mpi.h>]])
         # AC_COMPUTE_INT([OMPI_RELEASE], [OMPI_RELEASE_VERSION], [[#include <mpi.h>]])
         OMPI_MAJOR=`${GREP} OMPI_MAJOR_VERSION saved_conftest.i | cut -d' ' -f3`
         OMPI_MINOR=`${GREP} OMPI_MINOR_VERSION saved_conftest.i | cut -d' ' -f3`
         OMPI_RELEASE=`${GREP} OMPI_RELEASE_VERSION saved_conftest.i | cut -d' ' -f3`
         ompi_version="${OMPI_MAJOR}.${OMPI_MINOR}.${OMPI_RELEASE}"
         AC_MSG_RESULT(OpenMPI $ompi_version)
         unset OMPI_MAJOR
         unset OMPI_MINOR
         unset OMPI_RELEASE
         unset ompi_version
      fi
   else
      AC_MSG_RESULT([unknown])
   fi
   ${RM} -f saved_conftest.i
])

dnl It is not sufficient to use AC_CHECK_DECLS to check whether deprecated MPI
dnl constants are defined. OpenMPI still define the constants deprecated by MPI
dnl 3.0 as some error messages. We need to call AC_COMPILE_IFELSE to check
dnl whether they can be used without compilation errors.
dnl
m4_define([_UD_CHECK_MPI_CONSTANTS],
   [AC_MSG_CHECKING([whether $1 is defined])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([$4],[[int dummy=$1;]])],
                      [ac_have_const=yes], [ac_have_const=no])
    AC_MSG_RESULT([$ac_have_const])
    if test "x$ac_have_const" = xyes ; then
       AC_DEFINE([HAVE_$1], [1], [Define if $1 is defined and not deprecated])
    fi]
   [m4_ifvaln([$2$3], [AS_IF([test x$ac_have_const = xyes], [$2], [$3])])])

AC_DEFUN([UD_CHECK_MPI_CONSTANTS],
   [m4_map_args_sep([_$0(], [, [$2], [$3], [$4])], [], $1)])

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

# MPI_COMPILER_BASE
# -----------------------------------------------------------------
# Returns the base compiler command used in MPI compiler wrapper in
# ac_cv_mpi_compiler_base_MPICC, ac_cv_mpi_compiler_base_MPICXX,
# ac_cv_mpi_compiler_base_MPIF90, such as /usr/bin/gcc, /usr/bin/g++,
# /usr/bin/gfortran
#
AC_DEFUN([MPI_COMPILER_BASE],[
   # before checking, remove compile command-line options, if there is any
   compile_cmd=`echo $$1 | cut -d" " -f1`
   AC_MSG_CHECKING([base compiler command in $1 wrapper])
   compile_basename=
   case $compile_cmd in
        mpicc | mpicxx | mpif77 | mpif90 | mpifort | *[[\\/]]mpicc | *[[\\/]]mpicxx | *[[\\/]]mpif77 | *[[\\/]]mpif90 | *[[\\/]]mpifort )
           # MPICH, OpenMPI
           compile_basename=`$compile_cmd -show | cut -d' ' -f1`
           ;;
        mpixlc | mpixlcxx | mpixlf77 | mpixlf90 | *[[\\/]]mpixlc | *[[\\/]]mpixlcxx | *[[\\/]]mpixlf77 | *[[\\/]]mpixlf90 )
           # IBM XL MPI compilers
           compile_basename=`$compile_cmd -show | cut -d' ' -f1`
           ;;
        mpifccpx | mpiFCCpx | mpifrtpx | *[[\\/]]mpifccpx | *[[\\/]]mpiFCCpx | *[[\\/]]mpifrtpx )
           # Fujitsu MPI compilers: fccpx, FCCpx, frtpx
           compile_basename=`$compile_cmd -showme | cut -d' ' -f1`
           ;;
        cc | CC | ftn | *[[\\/]]cc | *[[\\/]]CC | *[[\\/]]ftn )
           # For Cray PrgEnv-intel, cc is a wrapper of icc
           # For Cray PrgEnv-gnu, cc is a wrapper of gcc
           eval "$compile_cmd --version" < /dev/null >& conftest.ver
           compile_basename=`head -n1 conftest.ver |cut -d' ' -f1`
           ${RM} -f conftest.ver
           if test "x${compile_basename}" = x ; then
              # For Cray PrgEnv-cray, cc is a wrapper of Cray CC
              # Cray cc -V sends the output to stderr.
              eval "$compile_cmd -V" < /dev/null >& conftest.ver
              compile_basename=`head -n1 conftest.ver |cut -d' ' -f1`
              ${RM} -f conftest.ver
           fi
           ;;
        *) break;;
   esac
   if test "x${compile_basename}" != x ; then
      AC_MSG_RESULT([$compile_basename])
      AC_PATH_PROG([ac_cv_mpi_compiler_base_$1], [$compile_basename])
   else
      AC_MSG_RESULT([not found])
   fi
   unset compile_basename
   unset compile_cmd
])# MPI_COMPILER_BASE

