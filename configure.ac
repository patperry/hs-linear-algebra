AC_PREREQ([2.61])
AC_INIT([Haskell blas], [0.9], [patperry@gmail.com])

AC_CONFIG_SRCDIR([blas.cabal])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH(compiler, 
    [AC_HELP_STRING([--with-compiler=<cmd>], [use Haskell compiler <cmd>])])

# Checks for programs
AC_PROG_CC_C99

# Checks for libraries
AC_CHECK_LIB([m],[cos])
AX_CHECK_BLAS

if test "$ax_blas_ok" = yes; then
  BUILD_PACKAGE_BOOL=True
else 
  BUILD_PACKAGE_BOOL=False
fi
AC_SUBST([BUILD_PACKAGE_BOOL])

AC_CONFIG_HEADERS([cbits/config.h])
AC_CONFIG_FILES([blas.buildinfo tests/Makefile])
AC_OUTPUT
