AC_PREREQ([2.61])
AC_INIT([Haskell Linear Algebra], [0.8], [patperry@gmail.com])

AC_CONFIG_SRCDIR([linear-algebra.cabal])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH(compiler, 
    [AC_HELP_STRING([--with-compiler=<cmd>], [use Haskell compiler <cmd>])])

# Checks for programs

# Checks for libraries
AC_CHECK_LIB([m],[cos])
AX_CHECK_LAPACK

if test "$ax_lapack_ok" = yes; then
  BUILD_PACKAGE_BOOL=True
else 
  BUILD_PACKAGE_BOOL=False
fi
AC_SUBST([BUILD_PACKAGE_BOOL])

AC_CONFIG_HEADERS([cbits/config.h])
AC_CONFIG_FILES([linear-algebra.buildinfo])
AC_OUTPUT
