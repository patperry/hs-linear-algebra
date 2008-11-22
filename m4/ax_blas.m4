
AC_DEFUN([AX_BLAS], [
AC_PREREQ(2.50)
ax_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

AX_BLAS_FUNC(gemm_, [ax_blas_ok=yes; ax_blas_underscore=yes],
    [AX_BLAS_FUNC(gemm, [ax_blas_ok=yes; ax_blas_underscore=no])])

AC_SUBST(BLAS_LIBS)

if test x"$ax_blas_ok" = xyes; then
    ifelse([$1],, [
        AC_DEFINE(HAVE_BLAS,1, [Define if you have a BLAS library.])
        AH_TEMPLATE([BLAS_FUNC], 
                [Define to a macro mangling the given BLAS function name])
        if test x"$ax_blas_underscore" = xyes; then
            AC_DEFINE([BLAS_FUNC(name)], [name ## _])
        else
            AC_DEFINE([BLAS_FUNC(name)], [name])
        fi
        ],
        [$1])
    :
else
    ax_blas_ok=no
    $2
fi

])dnl AX_BLAS
