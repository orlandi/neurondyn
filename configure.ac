AC_PREREQ(2.61)
AC_INIT([neruondyn], [0.1], [javierorlandi@javierorlandi.com],
        [neurondyntar], [http://www.javierorlandi.com/])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(config)
# Set this file to one that is in your directory
AC_CONFIG_SRCDIR([src/main.cc])
LT_INIT

# Automake initialization
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

# To enable CUDA
AC_ARG_ENABLE([cuda],
  AS_HELP_STRING([--enable-cuda],
    [Enable CUDA (default=yes)]),
    [enable_cuda=$enableval],
    [enable_cuda=yes])
if test "x$enable_cuda" = xyes; then
  AC_DEFINE(CUDA_ENABLED, 1, [Enable CUDA])
fi    
AM_CONDITIONAL(CUDA_ENABLED,
  [test "x$enable_cuda" = xyes])
# To enable OPENMP
AC_ARG_ENABLE([openmp],
  AS_HELP_STRING([--enable-openmp],
    [Enable OPENMP (default=yes)]),
    [enable_openmp=$enableval],
    [enable_openmp=yes])
if test "x$enable_openmp" = xyes; then
  AC_DEFINE(OPENMP_ENABLED, 1, [Enable OPENMP])
fi    
AM_CONDITIONAL(OPENMP_ENABLED,
  [test "x$enable_openmp" = xyes])

# Checks for programs
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
