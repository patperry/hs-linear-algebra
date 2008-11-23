
AC_PREREQ(2.61)
AC_INIT([Haskell blas], [0.7], [patperry@stanford.edu])
AC_CONFIG_SRCDIR([blas.cabal])

AC_CONFIG_HEADERS([cbits/config.h])

AC_PROG_CC
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH(compiler, 
    [AC_HELP_STRING([--with-compiler=<cmd>], [use Haskell compiler <cmd>])])

AX_BLAS
if test "$ax_blas_ok" = yes; then
  BUILD_PACKAGE_BOOL=True
else 
  BUILD_PACKAGE_BOOL=False
fi
AC_SUBST([BUILD_PACKAGE_BOOL])

AC_CONFIG_FILES([blas.buildinfo])
AC_OUTPUT
