dnl R interface to DSDP semidefinite programming library
dnl 
dnl By default the included DSDP code is built and linked, otherwise the user can supply
dnl the location of the dsdp library and includes, or the location of the dsdp installation directory
dnl
dnl 1) If the user sets PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS (all three) then these are used
dnl
dnl 2) If user provides location of the dsdp installation we use the lib and include directories there
dnl     --with-dsdp=<location of dsdp DIR>
dnl 3) If none of these are given the included DSDP code is built and linked
dnl
dnl Compiler flags for dsdp can be given using --with-dsdp-cflags

AC_INIT([Rdsdp], 0.1.1)

FOUND_ALL=""

dnl Check if global variables are already set, if so we are done
if test [-n "$PKG_CFLAGS" -a -n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS"] ; then
    FOUND_ALL="yes"
else
    if test [-n "$PKG_CFLAGS" -o -n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS"] ; then
	AC_MSG_ERROR([All three of PKG_CFLAGS PKG_CPPFLAGS and PKG_LIBS must be given together])
    fi
fi

dnl get arguments
AC_ARG_WITH([dsdp-cflags],
    [AS_HELP_STRING([--with-dsdp-cflags],[compiler flags for dsdp])],
    DSDP_CFLAGS=$withval,
    DSDP_CFLAGS="")

AC_ARG_WITH([dsdp-dir],
    [AS_HELP_STRING([--with-dsdp-dir],[specify location of dsdp installation])],
    DSDP_DIR=$withval,
    DSDP_DIR="")

if test [-z "$FOUND_ALL" -a -n "$DSDP_DIR"] ; then
    if test [! -d "$DSDP_DIR"] ; then
	AC_MSG_ERROR([DSDP installation directory $DSDP_DIR does not exist])
    fi

    PKG_CPPFLAGS="-I$DSDP_DIR/include"
    PKG_LIBS="-L$DSDP_DIR/lib"
    FOUND_ALL="yes"
fi

dnl No location given so build and link provided dsdp
if test [-z "$FOUND_ALL"] ; then
    echo "Using bundled DSDP"
    PKG_CPPFLAGS="-Idsdp/include"
    PKG_LIBS="-Ldsdp/lib"
    DSDP_TS="dsdp.ts"
else
    DSDP_TS=""
fi

dnl Find the R home directory.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "Could not determine R_HOME."
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

AC_PROG_CC
AC_PROG_CPP

PKG_LIBS="${PKG_LIBS} -ldsdp ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} -lm"
PKG_CFLAGS="-DNOSHORTS ${DSDP_CFLAGS}"

dnl if given by user, check header and library
if test [-n "$FOUND_ALL"] ; then
    LIBS="${PKG_LIBS}"
    CFLAGS="${CFLAGS} ${PKG_CFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"
    
    AC_CHECK_HEADER(declarations.h)
    if test "${ac_cv_header_declarations_h}" == no; then
	AC_MSG_ERROR([DSDP header declarations.h not found in provided location])
    fi
    
    AC_SEARCH_LIBS(easy_sdp,sdp,,AC_MSG_ERROR([Could not link DSDP library using $PKG_LIBS]))
fi

dnl substitute in src/Makevars.in
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(DSDP_TS)
AC_OUTPUT(src/Makevars)

exit 0

