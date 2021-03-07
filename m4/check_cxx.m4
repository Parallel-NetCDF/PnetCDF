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

